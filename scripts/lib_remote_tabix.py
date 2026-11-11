"""Read genomic windows from a remote BGZF/Tabix TSV using only stdlib + curl.

Coordinates passed to ``fetch`` are 1-based and inclusive. Returned strings retain
all columns and have no line terminator. No p-value filtering is performed.
Only generic TSV Tabix indexes are supported; SAM/VCF and CSI fail explicitly.
The cache contains original compressed bytes, response headers, checksums and an
index/header provenance manifest. ``offline=True`` forbids network access.

Specification: https://samtools.github.io/hts-specs/tabix.pdf
BGZF: https://samtools.github.io/hts-specs/SAMv1.pdf
"""

from __future__ import annotations

import hashlib
import concurrent.futures
import json
import re
import struct
import subprocess
import time
import threading
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit


# Shared by all source objects in this process, including each source's prefetch
# workers. Prevent several concurrently acquired tissues from flooding one host.
_NETWORK_LIMIT = threading.BoundedSemaphore(4)


class RemoteTabixError(RuntimeError):
    """Invalid remote response, damaged file, unsupported format, or cache miss."""


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def reg2bins(begin: int, end: int) -> list[int]:
    """Return TBI bins overlapping a 0-based, half-open interval."""
    if not 0 <= begin < end <= 1 << 29:
        raise ValueError("TBI interval must be nonempty and within [0, 2^29)")
    end -= 1
    bins = [0]
    for offset, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        bins.extend(range(offset + (begin >> shift), offset + (end >> shift) + 1))
    return bins


@dataclass(frozen=True)
class ReferenceIndex:
    bins: dict[int, tuple[tuple[int, int], ...]]
    linear: tuple[int, ...]


class TabixIndex:
    """Strict, little-endian TBI parser; virtual offsets stay Python integers."""

    def __init__(self, data: bytes):
        pos = 0

        def take(size: int) -> bytes:
            nonlocal pos
            if size < 0 or pos + size > len(data):
                raise RemoteTabixError("Truncated or invalid TBI index")
            result = data[pos:pos + size]
            pos += size
            return result

        def unpack(fmt: str):
            return struct.unpack("<" + fmt, take(struct.calcsize("<" + fmt)))

        def count() -> int:
            value = unpack("i")[0]
            if value < 0 or value > len(data):
                raise RemoteTabixError("Invalid count in TBI index")
            return value

        if take(4) != b"TBI\x01":
            raise RemoteTabixError("Expected TBI index (CSI is not supported)")
        nref = count()
        self.format, self.col_seq, self.col_beg, self.col_end, self.meta, self.skip = unpack("6i")
        names = take(count())
        if names and not names.endswith(b"\0"):
            raise RemoteTabixError("TBI reference names are not NUL terminated")
        self.names = tuple(x.decode("utf-8") for x in names[:-1].split(b"\0")) if names else ()
        if len(self.names) != nref or len(set(self.names)) != nref:
            raise RemoteTabixError("TBI reference count/name mismatch")
        if self.format not in (0, 0x10000):
            raise RemoteTabixError("Only generic TSV TBI indexes are supported")
        if self.col_seq < 1 or self.col_beg < 1 or self.col_end < 0 or self.skip < 0 or not 0 <= self.meta <= 255:
            raise RemoteTabixError("Invalid TBI column/header configuration")
        references = {}
        for name in self.names:
            bins = {}
            for _ in range(count()):
                bin_id = unpack("I")[0]
                nchunks = count()
                chunks = tuple(unpack("QQ") for _ in range(nchunks))
                if bin_id in bins:
                    raise RemoteTabixError("Duplicate TBI bin")
                # 37450 is a metadata pseudo-bin, whose pairs are not all offsets.
                if bin_id <= 37449:
                    if any(a >= b for a, b in chunks):
                        raise RemoteTabixError("Invalid TBI chunk bounds")
                    bins[bin_id] = chunks
            linear = unpack(str(count()) + "Q")
            references[name] = ReferenceIndex(bins, linear)
        if len(data) - pos not in (0, 8):
            raise RemoteTabixError("Unexpected trailing TBI bytes")
        self.references = references
        self.zero_based = bool(self.format & 0x10000)

    def chunks(self, chrom: str, begin: int, end: int) -> list[tuple[int, int]]:
        ref = self.references[chrom]
        interval = begin >> 14
        minimum = ref.linear[min(interval, len(ref.linear) - 1)] if ref.linear else 0
        chunks = sorted((max(a, minimum), b) for bin_id in reg2bins(begin, end)
                        for a, b in ref.bins.get(bin_id, ()) if b > minimum)
        merged: list[tuple[int, int]] = []
        for a, b in chunks:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        return merged


def _bgzf_size(header: bytes) -> int:
    if len(header) < 12 or header[:4] != b"\x1f\x8b\x08\x04":
        raise RemoteTabixError("Expected BGZF gzip member with FEXTRA flag")
    xlen = struct.unpack_from("<H", header, 10)[0]
    if xlen < 6 or len(header) < 12 + xlen:
        raise RemoteTabixError("Truncated BGZF extra header")
    pos, size = 12, None
    while pos < 12 + xlen:
        if pos + 4 > 12 + xlen:
            raise RemoteTabixError("Malformed BGZF extra subfield")
        tag, length = header[pos:pos + 2], struct.unpack_from("<H", header, pos + 2)[0]
        pos += 4
        if pos + length > 12 + xlen:
            raise RemoteTabixError("Malformed BGZF extra subfield length")
        if tag == b"BC":
            if length != 2 or size is not None:
                raise RemoteTabixError("Malformed or duplicate BGZF BC subfield")
            size = struct.unpack_from("<H", header, pos)[0] + 1
        pos += length
    if size is None or not 12 + xlen + 8 <= size <= 65536:
        raise RemoteTabixError("Invalid BGZF block size")
    return size


def bgzf_blocks(data: bytes, base_offset: int = 0):
    """Yield (compressed file offset, decompressed member), validating CRC/ISIZE."""
    offset = 0
    while offset < len(data):
        size = _bgzf_size(data[offset:])
        block = data[offset:offset + size]
        if len(block) != size:
            raise RemoteTabixError("Truncated BGZF member")
        if struct.unpack_from("<I", block, size - 4)[0] > 65536:
            raise RemoteTabixError("BGZF decompressed block exceeds 65536 bytes")
        try:
            decoder = zlib.decompressobj(wbits=31)  # Verify gzip CRC32 and ISIZE.
            raw = decoder.decompress(block, 65537)
        except (OSError, EOFError, ValueError, zlib.error) as exc:
            raise RemoteTabixError("BGZF CRC, length, or deflate validation failed") from exc
        if len(raw) > 65536:
            raise RemoteTabixError("BGZF decompressed block exceeds 65536 bytes")
        if not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
            raise RemoteTabixError("BGZF member boundary disagrees with its BSIZE")
        yield base_offset + offset, raw
        offset += size


def _response_fields(headers: bytes) -> tuple[bytes, dict[bytes, bytes]]:
    blocks = re.split(rb"\r?\n\r?\n", headers)
    responses = [part for part in blocks if re.match(rb"HTTP/\S+\s+\d{3}", part)]
    if not responses:
        raise RemoteTabixError("Missing HTTP response status")
    response = responses[-1]
    status = re.match(rb"HTTP/\S+\s+(\d{3})", response).group(1)
    fields = {}
    for line in response.splitlines()[1:]:
        if b":" in line:
            key, value = line.split(b":", 1)
            fields[key.strip().lower()] = value.strip()
    return status, fields


def validate_range_response(headers: bytes, body: bytes, start: int, end: int) -> tuple[int, int]:
    """Validate the final HTTP response, including redirect/proxy header chains."""
    status, fields = _response_fields(headers)
    if status != b"206":
        raise RemoteTabixError(f"Range request requires HTTP 206, received {status.decode()}; no full-download fallback")
    match = re.fullmatch(rb"bytes (\d+)-(\d+)/(\d+)", fields.get(b"content-range", b""))
    if not match:
        raise RemoteTabixError("Missing or invalid Content-Range")
    actual_start, actual_end, total = map(int, match.groups())
    expected_end = min(end, total - 1)
    if total < 1 or actual_start != start or actual_end != expected_end or len(body) != actual_end - actual_start + 1:
        raise RemoteTabixError("Content-Range or body length does not match requested bytes")
    if fields.get(b"content-encoding", b"identity").lower() not in (b"identity", b""):
        raise RemoteTabixError("Content-Encoding changed the BGZF byte offsets")
    length = fields.get(b"content-length")
    if length is not None:
        if not length.isdigit():
            raise RemoteTabixError("Invalid Content-Length")
        if (length.lstrip(b"0") or b"0") != str(len(body)).encode("ascii"):
            raise RemoteTabixError("Content-Length does not match downloaded range")
    return actual_end, total


class RemoteTabix:
    def __init__(self, url: str, cache_dir, *, index_url: str | None = None,
                 offline: bool = False, curl: str = "curl.exe", retries: int = 3,
                 max_span_bytes: int = 256 * 1024, max_index_bytes: int = 128 * 1024 * 1024):
        if urlsplit(url).scheme not in ("https", "http"):
            raise ValueError("RemoteTabix requires an HTTP(S) URL")
        if retries < 1 or max_span_bytes < 65536:
            raise ValueError("retries must be positive and max_span_bytes at least 65536")
        self.url = url
        self.index_url = index_url or url + ".tbi"
        self.offline, self.curl, self.retries = offline, curl, retries
        self.max_span_bytes = max_span_bytes
        self._lock = threading.RLock()
        self.cache_dir = Path(cache_dir) / _sha(url.encode())[:20]
        self.manifest_path = self.cache_dir / "manifest.json"
        if self.manifest_path.exists():
            self.manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            if self.manifest.get("url") != url or self.manifest.get("index_url") != self.index_url:
                raise RemoteTabixError("Cache belongs to a different source/index URL")
        else:
            if offline:
                raise RemoteTabixError("No cached manifest available for offline replay")
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.manifest = {"version": 1, "url": url, "index_url": self.index_url,
                             "ranges": [], "queries": [], "totals": {}}
        self._verified: dict[str, bytes] = {}
        self._header: list[str] | None = None
        self._index: TabixIndex | None = None
        self.max_index_bytes = max_index_bytes

    def _ensure_index(self):
        if self._index is not None:
            return
        first = self._read_range(self.index_url, 0, 65535)
        size = self.manifest["totals"][self.index_url]
        if size > self.max_index_bytes:
            raise RemoteTabixError(f"Index is {size} bytes, exceeding configured maximum")
        packed = first + (self._read_range(self.index_url, len(first), size - 1) if len(first) < size else b"")
        self._index = TabixIndex(b"".join(raw for _, raw in bgzf_blocks(packed)))
        if not self.offline:
            (self.cache_dir / "index.tbi").write_bytes(packed)
            self.manifest["index"] = {"file": "index.tbi", "sha256": _sha(packed),
                                      "size": len(packed), "references": list(self.references),
                                      "format": self.index.format, "col_seq": self.index.col_seq,
                                      "col_beg": self.index.col_beg, "col_end": self.index.col_end,
                                      "meta": self.index.meta, "skip": self.index.skip}
            self._save()

    @property
    def index(self) -> TabixIndex:
        self._ensure_index()
        return self._index

    @property
    def references(self) -> tuple[str, ...]:
        return self.index.names

    def _save(self):
        if not self.offline:
            with self._lock:
                temp = self.manifest_path.with_suffix(".json.tmp")
                temp.write_text(json.dumps(self.manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
                temp.replace(self.manifest_path)

    def _request(self, url: str, start: int, end: int) -> tuple[bytes, bytes]:
        """A bounded GET. Override in tests; no shell interpolation is used."""
        key = _sha(f"{url}\n{start}\n{end}".encode())
        output = self.cache_dir / (key + ".partial")
        headers = self.cache_dir / (key + ".headers.partial")
        command = [self.curl, "--ipv4", "--silent", "--show-error", "--location", "--fail",
                   "--connect-timeout", "20", "--max-time", "60",
                   "--speed-limit", "1024", "--speed-time", "20", "--max-filesize", str(end - start + 1),
                   "--header", "Accept-Encoding: identity", "--range", f"{start}-{end}",
                   "--dump-header", str(headers), "--output", str(output), url]
        try:
            with _NETWORK_LIMIT:
                result = subprocess.run(command, capture_output=True, check=False, timeout=75)
            if result.returncode:
                raise RemoteTabixError(f"curl range request failed ({result.returncode}): " + result.stderr.decode("utf-8", "replace")[-1500:])
            return headers.read_bytes(), output.read_bytes()
        finally:
            output.unlink(missing_ok=True)
            headers.unlink(missing_ok=True)

    def _download(self, url: str, start: int, end: int) -> dict:
        if self.offline:
            raise RemoteTabixError(f"Offline cache miss for {url} bytes {start}-{end}")
        last_error = None
        for attempt in range(self.retries):
            try:
                headers, body = self._request(url, start, end)
                actual_end, total = validate_range_response(headers, body, start, end)
                key = _sha(f"{url}\n{start}\n{end}".encode())
                filename = key + ".bin"
                entry = {"url": url, "requested_start": start, "requested_end": end,
                         "start": start, "end": actual_end, "total": total,
                         "size": len(body), "sha256": _sha(body), "file": filename,
                         "headers": key + ".headers", "headers_sha256": _sha(headers),
                         "downloaded_utc": datetime.now(timezone.utc).isoformat()}
                with self._lock:
                    self._check_validators(url, headers)
                    known = self.manifest["totals"].get(url)
                    if known is not None and known != total:
                        raise RemoteTabixError("Remote file size changed within cached acquisition")
                    for name, content in ((filename, body), (key + ".headers", headers)):
                        temporary = self.cache_dir / (name + ".tmp")
                        temporary.write_bytes(content)
                        temporary.replace(self.cache_dir / name)
                    self.manifest["ranges"] = [x for x in self.manifest["ranges"] if x['file'] != filename]
                    self.manifest["ranges"].append(entry)
                    self.manifest["totals"][url] = total
                    self._save()
                return entry
            except (RemoteTabixError, OSError, subprocess.TimeoutExpired) as exc:
                last_error = exc
                if attempt + 1 < self.retries:
                    time.sleep(min(2 ** attempt, 8))
        raise RemoteTabixError(f"Failed after {self.retries} range attempts: {last_error}") from last_error

    def _check_validators(self, url: str, headers: bytes):
        _, fields = _response_fields(headers)
        values = {key.decode(): fields[key].decode("latin-1")
                  for key in (b"etag", b"last-modified") if key in fields}
        known = self.manifest.setdefault("validators", {}).get(url)
        if known is not None:
            for key, value in known.items():
                if values.get(key) != value:
                    raise RemoteTabixError(f"Remote {key} changed or disappeared within cached acquisition")
        else:
            self.manifest["validators"][url] = values

    def _cached_bytes(self, entry: dict) -> bytes:
        name = entry["file"]
        if name not in self._verified:
            data = (self.cache_dir / name).read_bytes()
            if len(data) != entry["size"] or _sha(data) != entry["sha256"]:
                raise RemoteTabixError(f"Cached range checksum mismatch: {name}")
            headers = (self.cache_dir / entry["headers"]).read_bytes()
            if entry.get("headers_sha256") and _sha(headers) != entry["headers_sha256"]:
                raise RemoteTabixError(f"Cached response header checksum mismatch: {name}")
            last, total = validate_range_response(headers, data, entry["requested_start"], entry["requested_end"])
            if last != entry["end"] or total != entry["total"] or total != self.manifest["totals"][entry["url"]]:
                raise RemoteTabixError("Cached range metadata disagrees with its HTTP response")
            self._check_validators(entry["url"], headers)
            self._verified[name] = data
            # Keep decompressed extraction streaming; compressed memory is bounded.
            while sum(map(len, self._verified.values())) > 2 * self.max_span_bytes and len(self._verified) > 1:
                del self._verified[next(iter(self._verified))]
        return self._verified[name]

    def _read_range(self, url: str, start: int, end: int) -> bytes:
        if start < 0 or end < start:
            raise ValueError("Invalid compressed byte interval")
        known = self.manifest["totals"].get(url)
        if known is None:
            # Establish the exact object size before parallel requests; this also
            # handles a short index whose initial 64 KiB range is clipped at EOF.
            self._download(url, start, min(end, start + self.max_span_bytes - 1))
            known = self.manifest["totals"][url]
        end = min(end, known - 1)
        entries = sorted((x for x in self.manifest["ranges"] if x["url"] == url), key=lambda x: x["start"])
        gaps = []
        cursor = start
        for entry in entries:
            if entry["end"] < cursor:
                continue
            if entry["start"] > end:
                break
            while cursor < entry["start"]:
                gap_end = min(end, entry["start"] - 1, cursor + self.max_span_bytes - 1)
                gaps.append((cursor, gap_end))
                cursor = gap_end + 1
            cursor = max(cursor, entry["end"] + 1)
        while cursor <= end:
            gap_end = min(end, cursor + self.max_span_bytes - 1)
            gaps.append((cursor, gap_end))
            cursor = gap_end + 1
        if self.offline and gaps:
            raise RemoteTabixError(f"Offline cache miss for {url} bytes {gaps[0][0]}-{gaps[0][1]}")
        # Keep at most four pending requests per object and four in flight across
        # all objects. A failed batch never starts later batches; successful
        # sibling ranges remain available to the next resumable run.
        if len(gaps) == 1:
            self._download(url, *gaps[0])
        elif gaps:
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                for batch_start in range(0, len(gaps), 4):
                    jobs = [pool.submit(self._download, url, a, b) for a, b in gaps[batch_start:batch_start + 4]]
                    for job in jobs:
                        job.result()
        pieces = []
        pos = start
        while pos <= end:
            entries = [x for x in self.manifest["ranges"] if x["url"] == url]
            covering = [x for x in entries if x["start"] <= pos <= x["end"]]
            if covering:
                entry = max(covering, key=lambda x: x["end"])
            else:
                raise RemoteTabixError("Downloaded range cache does not cover requested bytes")
            last = min(end, entry["end"])
            pieces.append(self._cached_bytes(entry)[pos - entry["start"]:last - entry["start"] + 1])
            pos = last + 1
        return b"".join(pieces)

    def _block_size(self, offset: int) -> int:
        head = self._read_range(self.url, offset, offset + 17)
        if len(head) < 12:
            raise RemoteTabixError("BGZF offset lies beyond the remote file")
        xlen = struct.unpack_from("<H", head, 10)[0]
        if xlen > 65512:
            raise RemoteTabixError("Invalid BGZF extra header length")
        if 12 + xlen > len(head):
            head = self._read_range(self.url, offset, offset + 11 + xlen)
        return _bgzf_size(head)

    def get_header_lines(self) -> list[str]:
        """Return initial skipped/comment/header lines, without newline characters."""
        if self._header is not None:
            return list(self._header)
        result, carry, offset, line_number = [], b"", 0, 0
        consumed = 0
        while consumed < 4 * 1024 * 1024:
            size = self._block_size(offset)
            packed = self._read_range(self.url, offset, offset + size - 1)
            raw = next(bgzf_blocks(packed))[1]
            consumed += len(raw)
            lines = (carry + raw).split(b"\n")
            carry = lines.pop()
            if not raw and carry:
                lines.append(carry)
                carry = b""
            for line in lines:
                text = line.rstrip(b"\r").decode("utf-8")
                line_number += 1
                is_header = line_number <= self.index.skip or self._comment(text)
                if not is_header and line_number == 1:
                    try:
                        self._record(text)
                    except (RemoteTabixError, ValueError):
                        is_header = True
                if is_header:
                    result.append(text)
                elif text:
                    self._header = result
                    if not self.offline:
                        self.manifest["header"] = {"lines": result, "compressed_prefix_end": offset + size - 1}
                        self._save()
                    return list(result)
            if not raw:
                self._header = result
                return list(result)
            offset += size
        raise RemoteTabixError("Header exceeds the configured 4 MiB decompressed limit")

    @property
    def header(self) -> str:
        """Last header line, usually the TSV column names (empty if absent)."""
        return self.get_header()

    def get_header(self) -> str:
        """Return the TSV column header; all header lines via get_header_lines()."""
        lines = self.get_header_lines()
        return lines[-1] if lines else ""

    def _comment(self, line: str) -> bool:
        return bool(self.index.meta and line.startswith(chr(self.index.meta)))

    def _record(self, line: str) -> tuple[str, int, int]:
        fields = line.split("\t")
        try:
            chrom = fields[self.index.col_seq - 1]
            start = int(fields[self.index.col_beg - 1])
            end = int(fields[self.index.col_end - 1]) if self.index.col_end else start
            if self.index.zero_based:
                if not self.index.col_end or self.index.col_end == self.index.col_beg:
                    end = start + 1
            else:
                start -= 1
            if start < 0 or end < start:
                raise ValueError("Invalid genomic coordinates")
            # htslib represents a zero-length BED feature as a one-base interval.
            end = max(end, start + 1)
            return chrom, start, end
        except (IndexError, ValueError) as exc:
            raise RemoteTabixError("Malformed indexed TSV record: " + line[:160]) from exc

    def _chunk_data(self, begin: int, end: int) -> bytes:
        first_block, first_inner = begin >> 16, begin & 65535
        last_block, last_inner = end >> 16, end & 65535
        stop = last_block + self._block_size(last_block) if last_inner else last_block
        if stop <= first_block:
            raise RemoteTabixError("Invalid chunk compressed boundaries")
        packed = self._read_range(self.url, first_block, stop - 1)
        parts = []
        for offset, raw in bgzf_blocks(packed, first_block):
            lo = first_inner if offset == first_block else 0
            hi = last_inner if offset == last_block else len(raw)
            if lo > len(raw) or hi > len(raw) or hi < lo:
                raise RemoteTabixError("Virtual offset exceeds BGZF uncompressed block")
            parts.append(raw[lo:hi])
        return b"".join(parts)

    def fetch(self, chrom: str, start: int, end: int):
        """Yield every indexed row overlapping a 1-based inclusive window."""
        if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start <= end <= 1 << 29:
            raise ValueError("fetch expects 1-based inclusive coordinates within the TBI limit")
        chrom = str(chrom)
        if chrom not in self.index.references:
            alias = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alias not in self.index.references:
                raise KeyError(f"Reference {chrom!r} is not present in the TBI index")
            chrom = alias
        begin = start - 1
        chunks = self.index.chunks(chrom, begin, end)
        header_lines = set(self.get_header_lines())
        count = 0
        for a, b in chunks:
            text = self._chunk_data(a, b).decode("utf-8")
            for line in text.splitlines():
                if not line or self._comment(line) or line in header_lines:
                    continue
                ref, record_start, record_end = self._record(line)
                if ref == chrom and record_start < end and record_end > begin:
                    count += 1
                    yield line
        if not self.offline:
            self.manifest["queries"].append({"chrom": chrom, "start": start, "end": end,
                                              "chunks": [[a, b] for a, b in chunks], "records": count})
            self._save()
