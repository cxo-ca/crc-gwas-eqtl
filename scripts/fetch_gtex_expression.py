#!/usr/bin/env python3
"""Fetch two GTEx v8 expression members by HTTP range and measure phenotype SD.

Only Python's standard library and curl are required. The source archive is an
uncompressed TAR containing compressed BED members. Every requested range has a
checked HTTP response, SHA-256, and JSON provenance record, permitting --offline
replay. This script does not download the full archive or impute phenotype SD.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import statistics
import subprocess
import tarfile
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/"
    "single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar"
)
TISSUES = ("Colon_Sigmoid", "Colon_Transverse")
SOURCE_NOTES = {
    "download_documentation": "https://gtexportal.org/home/downloads/adult-gtex/overview",
    "normalization_documentation": "https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md",
    "normalization_code": "https://github.com/broadinstitute/pyqtl/blob/master/qtl/norm.py",
    "interpretation": (
        "Sample standard deviation (denominator N-1) of the published fully "
        "processed, filtered, normalized expression values used as FastQTL input. "
        "TMM normalization precedes inverse-normal transformation across donors. "
        "These SDs are measured before FastQTL covariate residualization; they "
        "are not assumed to be 1 and are not calculated from raw TPM."
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def checked_content_range(raw_headers: str, start: int, end: int) -> tuple[int, dict[str, str]]:
    """Validate the last HTTP response, following redirects and proxy handshakes."""
    responses = re.split(r"\r?\n\r?\n", raw_headers.strip())
    responses = [item for item in responses if re.match(r"HTTP/\S+\s+\d+", item)]
    if not responses:
        raise ValueError("No HTTP response headers were recorded")
    lines = responses[-1].splitlines()
    status = int(lines[0].split()[1])
    if status != 206:
        raise ValueError(f"Expected HTTP 206 for a byte range, got {status}; refusing full archive")
    headers = {}
    for line in lines[1:]:
        if ":" in line:
            key, value = line.split(":", 1)
            headers[key.lower().strip()] = value.strip()
    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", headers.get("content-range", ""))
    if not match:
        raise ValueError("Missing or invalid Content-Range")
    actual_start, actual_end, total = map(int, match.groups())
    if (actual_start, actual_end) != (start, end) or total <= end:
        raise ValueError(f"Content-Range does not match requested bytes {start}-{end}")
    if headers.get("content-encoding", "identity").lower() != "identity":
        raise ValueError("Unexpected content encoding would invalidate byte offsets")
    if "content-length" in headers and int(headers["content-length"]) != end - start + 1:
        raise ValueError("Content-Length does not match requested range")
    return total, headers


class RangeCache:
    def __init__(self, url: str, directory: Path, offline: bool, curl: str):
        self.url = url
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.curl = curl
        self.total_size: int | None = None
        self.etag: str | None = None
        self.used: list[dict] = []

    def get(self, start: int, length: int, purpose: str) -> Path:
        if start < 0 or length <= 0:
            raise ValueError("Invalid requested byte range")
        end = start + length - 1
        prefix = f"{start:012d}-{end:012d}"
        body = self.directory / (prefix + ".bin")
        headers_file = self.directory / (prefix + ".headers.txt")
        record_file = self.directory / (prefix + ".json")
        if record_file.exists():
            record = json.loads(record_file.read_text(encoding="utf-8"))
            if record["url"] != self.url or record["start"] != start or record["end"] != end:
                raise ValueError(f"Cached provenance mismatch: {record_file}")
            if not body.exists() or body.stat().st_size != length or sha256_file(body) != record["sha256"]:
                raise ValueError(f"Corrupt cached byte range: {body}")
            if not headers_file.exists() or sha256_file(headers_file) != record["headers_sha256"]:
                raise ValueError(f"Missing or altered cached response headers: {headers_file}")
            total, headers = checked_content_range(headers_file.read_text(encoding="iso-8859-1"), start, end)
            if total != record["archive_size"]:
                raise ValueError("Cached source size mismatch")
        else:
            if self.offline:
                raise FileNotFoundError(f"Offline range missing: {body}")
            temporary = body.with_name(body.name + ".part")
            temporary_headers = headers_file.with_name(headers_file.name + ".part")
            command = [
                self.curl, "--fail", "--location", "--silent", "--show-error",
                "--retry", "3", "--retry-all-errors", "--connect-timeout", "20",
                "--max-time", "600", "--max-filesize", str(length),
                "--header", "Accept-Encoding: identity", "--range", f"{start}-{end}",
                "--dump-header", str(temporary_headers), "--output", str(temporary), self.url,
            ]
            print(f"Fetch {purpose}: bytes {start}-{end} ({length:,} bytes)", flush=True)
            subprocess.run(command, check=True)
            total, headers = checked_content_range(temporary_headers.read_text(encoding="iso-8859-1"), start, end)
            if temporary.stat().st_size != length:
                raise ValueError("Downloaded byte count differs from requested range")
            record = {
                "url": self.url, "start": start, "end": end,
                "bytes": length, "archive_size": total,
                "downloaded_at_utc": utc_now(), "purpose": purpose,
                "http_status": 206, "content_range": headers["content-range"],
                "etag": headers.get("etag"), "last_modified": headers.get("last-modified"),
                "sha256": sha256_file(temporary), "headers_sha256": sha256_file(temporary_headers),
                "data_file": str(body.resolve()), "headers_file": str(headers_file.resolve()),
            }
            temporary.replace(body)
            temporary_headers.replace(headers_file)
            write_json(record_file, record)
        if self.total_size is not None and self.total_size != total:
            raise ValueError("Source archive size changed between range requests")
        current_etag = headers.get("etag")
        if self.etag is not None and self.etag != current_etag:
            raise ValueError("Source archive ETag changed between range requests")
        self.total_size = total
        self.etag = current_etag
        self.used.append({**record, "record_file": str(record_file.resolve())})
        return body


def parse_pax(payload: bytes) -> dict[str, str]:
    """Parse POSIX length-prefixed PAX records (length counts bytes)."""
    result: dict[str, str] = {}
    offset = 0
    while offset < len(payload):
        space = payload.find(b" ", offset)
        if space < 0 or not payload[offset:space].isdigit():
            raise ValueError("Malformed PAX record length")
        size = int(payload[offset:space])
        finish = offset + size
        if size <= space - offset + 1 or finish > len(payload) or payload[finish - 1:finish] != b"\n":
            raise ValueError("Malformed PAX record bounds")
        content = payload[space + 1:finish - 1]
        if b"=" not in content:
            raise ValueError("Malformed PAX key/value")
        key, value = content.split(b"=", 1)
        result[key.decode("utf-8")] = value.decode("utf-8")
        offset = finish
    return result


def target_tissue(name: str) -> str | None:
    base = PurePosixPath(name).name
    for tissue in TISSUES:
        if re.fullmatch(re.escape(tissue) + r"(?:\.[^.]+)*\.(?:normalized_)?expression\.bed\.gz", base):
            return tissue
    return None


def fetch_members(cache: RangeCache, raw_directory: Path) -> list[dict]:
    offset = 0
    found = {}
    global_pax: dict[str, str] = {}
    pending_pax: dict[str, str] = {}
    long_name: str | None = None
    for _ in range(1000):
        header_path = cache.get(offset, 512, "tar_header")
        block = header_path.read_bytes()
        if block == b"\0" * 512:
            # Archive traversal stops early when both requested members are found.
            # If reached otherwise, require both TAR terminator blocks.
            second = cache.get(offset + 512, 512, "tar_terminator").read_bytes()
            if second != b"\0" * 512:
                raise ValueError("Invalid TAR end marker")
            break
        info = tarfile.TarInfo.frombuf(block, "utf-8", "strict")  # validates header checksum
        payload_offset = offset + 512
        if info.size < 0:
            raise ValueError("Negative TAR member size")
        if info.type in (tarfile.XHDTYPE, tarfile.XGLTYPE, tarfile.GNUTYPE_LONGNAME, tarfile.GNUTYPE_LONGLINK):
            if info.size > 1024 * 1024:
                raise ValueError("Unexpectedly large TAR extension")
            payload = cache.get(payload_offset, info.size, "tar_extension").read_bytes() if info.size else b""
            if info.type == tarfile.XHDTYPE:
                pending_pax.update(parse_pax(payload))
            elif info.type == tarfile.XGLTYPE:
                global_pax.update(parse_pax(payload))
            elif info.type == tarfile.GNUTYPE_LONGNAME:
                long_name = payload.rstrip(b"\0").decode("utf-8")
            offset = payload_offset + ((info.size + 511) // 512) * 512
            continue
        metadata = {**global_pax, **pending_pax}
        name = metadata.get("path", long_name or info.name)
        size = int(metadata.get("size", info.size))
        if size < 0 or (cache.total_size is not None and payload_offset + size > cache.total_size):
            raise ValueError("TAR member extends beyond source archive")
        tissue = target_tissue(name)
        if tissue:
            if not info.isreg() or any(key.startswith("GNU.sparse") for key in metadata):
                raise ValueError(f"Target must be a regular non-sparse member: {name}")
            if tissue in found:
                raise ValueError(f"Duplicate target tissue member: {name}")
            if not 0 < size <= 250 * 1024 * 1024:
                raise ValueError(f"Unexpected expression member size: {size}")
            payload_path = cache.get(payload_offset, size, f"expression:{tissue}")
            raw_path = raw_directory / PurePosixPath(name).name
            if raw_path.exists() and sha256_file(raw_path) != sha256_file(payload_path):
                raise ValueError(f"Existing raw member has different content: {raw_path}")
            if not raw_path.exists():
                shutil.copyfile(payload_path, raw_path)
            found[tissue] = {
                "tissue": tissue, "member_name": name, "header_offset": offset,
                "payload_offset": payload_offset, "bytes": size,
                "sha256": sha256_file(raw_path), "raw_file": str(raw_path.resolve()),
                "header_checksum_verified": True,
            }
            if len(found) == len(TISSUES):
                return [found[tissue] for tissue in TISSUES]
        offset = payload_offset + ((size + 511) // 512) * 512
        pending_pax = {}
        long_name = None
    missing = sorted(set(TISSUES) - set(found))
    raise ValueError(f"Requested expression BED members not found: {missing}")


def read_candidates(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise ValueError("Candidate file is missing gene_id column")
        result = {}
        for row in reader:
            gene_id = row["gene_id"].strip()
            if not re.fullmatch(r"ENSG\d+(?:\.\d+)?", gene_id):
                raise ValueError(f"Unexpected candidate gene ID: {gene_id}")
            stable = gene_id.split(".", 1)[0]
            if stable in result and result[stable] != gene_id:
                raise ValueError(f"Conflicting candidate gene versions: {stable}")
            result[stable] = gene_id
    if not result:
        raise ValueError("Candidate gene list is empty")
    return result


def extract_sd(members: list[dict], candidates: dict[str, str]) -> tuple[list[dict], list[dict]]:
    rows = []
    coverage = []
    for member in members:
        found = set()
        with gzip.open(member["raw_file"], "rt", encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = next(reader)
            if len(header) < 6 or header[0].lstrip("#").lower() not in ("chr", "chrom", "chromosome"):
                raise ValueError("Unexpected GTEx expression BED header")
            if header[1:3] != ["start", "end"] or header[3] not in ("gene_id", "phenotype_id"):
                raise ValueError("Unexpected GTEx expression BED coordinate/ID columns")
            n = len(header) - 4
            if len(set(header[4:])) != n:
                raise ValueError("Duplicate expression sample IDs")
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(header):
                    raise ValueError(f"Expression BED field count mismatch at line {line_number}")
                gene_id = row[3]
                stable = gene_id.split(".", 1)[0]
                if stable not in candidates:
                    continue
                if gene_id != candidates[stable]:
                    raise ValueError(f"Candidate/expression gene version mismatch: {candidates[stable]} vs {gene_id}")
                if stable in found:
                    raise ValueError(f"Duplicate candidate phenotype: {member['tissue']} {gene_id}")
                values = [float(value) for value in row[4:]]
                if not all(math.isfinite(value) for value in values):
                    raise ValueError(f"Non-finite phenotype value: {gene_id}")
                sd = statistics.stdev(values)
                start, end = int(row[1]), int(row[2])
                if not math.isfinite(sd) or sd <= 0 or start < 0 or end != start + 1:
                    raise ValueError(f"Invalid SD or non-point TSS BED interval: {gene_id}")
                rows.append({
                    "gene_id": gene_id, "candidate_gene_id": candidates[stable],
                    "gene_id_stable": stable, "tissue": member["tissue"],
                    "n": n, "sdY": format(sd, ".17g"), "chr": row[0],
                    "bed_start": start, "bed_end": end, "tss": end,
                    "genome_build": "GRCh38", "sdY_method": "sample_sd_of_published_normalized_expression",
                    "source_url": SOURCE_URL, "source_member": member["member_name"],
                    "source_member_sha256": member["sha256"],
                })
                found.add(stable)
        coverage.append({
            "tissue": member["tissue"], "n": n, "candidate_count": len(candidates),
            "observed_count": len(found),
            "missing_candidate_gene_ids": [candidates[stable] for stable in sorted(set(candidates) - found)],
        })
    rows.sort(key=lambda row: (row["tissue"], row["gene_id"]))
    return rows, coverage


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true", help="Use only previously validated range caches")
    parser.add_argument("--candidates", type=Path, default=ROOT / "results/public/candidate_genes.tsv")
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "data/raw/summary_statistics/GTEx_v8/expression")
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/acquisition/gtex_expression_sdY.tsv")
    parser.add_argument("--curl", default="curl.exe" if os.name == "nt" else "curl")
    args = parser.parse_args()
    candidates = read_candidates(args.candidates)
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    cache = RangeCache(SOURCE_URL, args.raw_dir / "ranges", args.offline, args.curl)
    members = fetch_members(cache, args.raw_dir)
    rows, coverage = extract_sd(members, candidates)
    if not rows:
        raise ValueError("No candidate expression phenotypes were found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = args.output.with_name(args.output.name + ".tmp")
    with temporary_output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary_output.replace(args.output)
    manifest = {
        "completed_at_utc": utc_now(), "offline_replay": args.offline,
        "source_url": SOURCE_URL, "source_archive_size": cache.total_size,
        "source_etag": cache.etag, "source_notes": SOURCE_NOTES,
        "traversal": "Validated TAR headers until both requested members were located; no full archive download",
        "bed_coordinates": "0-based half-open point intervals; tss = bed_end is 1-based GRCh38",
        "candidate_file": str(args.candidates.resolve()), "candidate_sha256": sha256_file(args.candidates),
        "output_file": str(args.output.resolve()), "output_sha256": sha256_file(args.output),
        "output_rows": len(rows), "members": members, "coverage": coverage,
        "range_count": len(cache.used), "selected_range_bytes": sum(item["bytes"] for item in cache.used),
        "ranges": cache.used,
    }
    write_json(args.raw_dir / "acquisition_manifest.json", manifest)
    print(json.dumps({"output": str(args.output), "rows": len(rows), "coverage": coverage}, indent=2), flush=True)


if __name__ == "__main__":
    main()
