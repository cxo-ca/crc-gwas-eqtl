"""Offline fixtures exercise real BGZF members and binary TBI virtual offsets."""

import json
import struct
import sys
import tempfile
import unittest
import zlib
import concurrent.futures
import threading
import time
from types import SimpleNamespace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lib_remote_tabix import (RemoteTabix, RemoteTabixError, TabixIndex,
                              bgzf_blocks, reg2bins, validate_range_response)


URL = "https://example.invalid/study.tsv.gz"


def bgzf(raw):
    compressor = zlib.compressobj(wbits=-15)
    compressed = compressor.compress(raw) + compressor.flush()
    size = 18 + len(compressed) + 8
    assert size <= 65536 and len(raw) <= 65536
    return (b"\x1f\x8b\x08\x04" + struct.pack("<I", 0) + b"\x00\xff"
            + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, size - 1)
            + compressed + struct.pack("<II", zlib.crc32(raw), len(raw)))


def smallest_bin(start, end):
    return reg2bins(start, end)[-1] if start >> 14 == (end - 1) >> 14 else next(
        offset + (start >> shift) for offset, shift in
        ((585, 17), (73, 20), (9, 23), (1, 26), (0, 29))
        if start >> shift == (end - 1) >> shift)


def fixture(*, bed=False):
    if bed:
        header = "#chrom\tstart\tend\tvalue"
        rows = ["chr1\t0\t1\ta", "chr1\t89\t110\tb", "chr1\t110\t111\tc",
                "chr1\t16380\t16390\td", "chr1\t40000\t40001\te"]
        parsed = [(r.split("\t")[0], int(r.split("\t")[1]), int(r.split("\t")[2])) for r in rows]
    else:
        header = "#chrom\tposition\tp_value"
        rows = ["chr1\t1\t0.8", "chr1\t2\t0.9", "chr1\t2\t0.9", "chr1\t16384\t0.5",
                "chr1\t16385\t0.2", "chr1\t32768\t0.3", "chr1\t40000\t" + "x" * 65000,
                "chr1\t50000\t0.7", "chr2\t2\t0.6"]
        parsed = [(r.split("\t")[0], int(r.split("\t")[1]) - 1, int(r.split("\t")[1])) for r in rows]
    raw = (header + "\n" + "\n".join(rows) + "\n").encode()
    # Deliberately split one line and end another exactly at a BGZF boundary.
    record_starts = []
    cursor = len(header.encode()) + 1
    for row in rows:
        record_starts.append(cursor)
        cursor += len(row.encode()) + 1
    cuts = sorted(set([0, record_starts[1], record_starts[3] + 4,
                       min(32000, len(raw)), min(64000, len(raw)), len(raw)]))
    blocks, locations, compressed_cursor = [], [], 0
    for left, right in zip(cuts, cuts[1:]):
        packed = bgzf(raw[left:right])
        locations.append((left, right, compressed_cursor))
        blocks.append(packed)
        compressed_cursor += len(packed)

    def virtual(pos):
        if pos == len(raw):
            return compressed_cursor << 16
        for left, right, offset in locations:
            if left <= pos < right:
                return (offset << 16) | (pos - left)
        raise AssertionError(pos)

    names = tuple(dict.fromkeys(chrom for chrom, _, _ in parsed))
    names_bytes = b"".join(name.encode() + b"\0" for name in names)
    tbi = b"TBI\x01" + struct.pack("<8i", len(names), 0x10000 if bed else 0,
                                   1, 2, 3 if bed else 2, ord("#"), 0, len(names_bytes)) + names_bytes
    for name in names:
        bins = {}
        selected = []
        for idx, (chrom, start, end) in enumerate(parsed):
            if chrom != name:
                continue
            a, b = virtual(record_starts[idx]), virtual(record_starts[idx] + len(rows[idx].encode()) + 1)
            bins.setdefault(smallest_bin(start, end), []).append((a, b))
            selected.append((start, end, a))
        # Duplicate overlapping index chunks must not duplicate rows on output.
        first_bin = next(iter(bins))
        bins[first_bin].append(bins[first_bin][0])
        tbi += struct.pack("<i", len(bins))
        for bin_id, chunks in bins.items():
            tbi += struct.pack("<Ii", bin_id, len(chunks))
            tbi += b"".join(struct.pack("<QQ", *chunk) for chunk in chunks)
        nintervals = 1 + (max(end for _, end, _ in selected) >> 14)
        linear = []
        previous = 0
        for interval in range(nintervals):
            overlapping = [offset for start, end, offset in selected
                           if start < (interval + 1) * 16384 and end > interval * 16384]
            previous = min(overlapping) if overlapping else previous
            linear.append(previous)
        tbi += struct.pack("<i", len(linear)) + struct.pack("<" + "Q" * len(linear), *linear)
    tbi += struct.pack("<Q", 0)
    return {URL: b"".join(blocks) + bgzf(b""), URL + ".tbi": bgzf(tbi) + bgzf(b"")}, header, rows, parsed


class MemoryRemote(RemoteTabix):
    def __init__(self, cache_dir, blobs, **kwargs):
        self.blobs = blobs
        self.requests = []
        self.etag = '"fixture-v1"'
        super().__init__(URL, cache_dir, retries=1, **kwargs)

    def _request(self, url, start, end):
        self.requests.append((url, start, end))
        blob = self.blobs[url]
        body = blob[start:end + 1]
        actual_end = start + len(body) - 1
        header = (f"HTTP/1.1 200 Connection established\r\n\r\n"
                  f"HTTP/1.1 302 Found\r\nLocation: {url}\r\n\r\n"
                  f"HTTP/2 206\r\nContent-Range: bytes {start}-{actual_end}/{len(blob)}\r\n"
                  f"ETag: {self.etag}\r\nLast-Modified: Tue, 22 Sep 2026 00:00:00 GMT\r\n"
                  f"Content-Length: {len(body)}\r\n\r\n").encode()
        return header, body


class RemoteTabixTests(unittest.TestCase):
    def test_lazy_header_and_point_boundaries(self):
        blobs, header, rows, parsed = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            self.assertEqual(source.requests, [])
            self.assertEqual(source.references, ("chr1", "chr2"))
            self.assertEqual(source.get_header(), header)
            self.assertEqual(source.get_header_lines(), [header])
            for start, end in [(1, 1), (2, 2), (3, 16383), (16384, 16384), (16385, 16385),
                               (16384, 16385), (32768, 40000), (50000, 50000), (50001, 60000), (1, 50000)]:
                expected = [row for row, (chrom, a, b) in zip(rows, parsed)
                            if chrom == "chr1" and a < end and b > start - 1]
                self.assertEqual(list(source.fetch("1", start, end)), expected, (start, end))
            self.assertEqual(list(source.fetch("chr2", 2, 2)), [rows[-1]])
            # Identical duplicate source rows are real records and are preserved.
            self.assertEqual(list(source.fetch("chr1", 2, 2)), rows[1:3])
            self.assertTrue(all(end - start + 1 <= source.max_span_bytes for _, start, end in source.requests))

    def test_bed_overlap_coordinates(self):
        blobs, _, rows, parsed = fixture(bed=True)
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            for start, end in [(1, 1), (89, 89), (90, 90), (100, 100), (110, 110), (111, 111),
                               (16384, 16385), (16390, 16390), (16391, 16391)]:
                expected = [row for row, (_, a, b) in zip(rows, parsed) if a < end and b > start - 1]
                self.assertEqual(list(source.fetch("chr1", start, end)), expected, (start, end))

    def test_offline_replay_checksums_and_cache_reuse(self):
        blobs, header, _, _ = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            expected = list(source.fetch("chr1", 1, 50000))
            nrequests = len(source.requests)
            self.assertEqual(list(source.fetch("chr1", 2, 40000)), expected[1:-1])
            self.assertEqual(len(source.requests), nrequests)
            offline = RemoteTabix(URL, directory, offline=True)
            self.assertEqual(offline.get_header(), header)
            self.assertEqual(list(offline.fetch("chr1", 1, 50000)), expected)
            manifest = json.loads(source.manifest_path.read_text())
            self.assertEqual(manifest["queries"][0]["records"], len(expected))
            self.assertTrue((source.cache_dir / manifest["index"]["file"]).exists())
            entry = manifest["ranges"][0]
            (source.cache_dir / entry["file"]).write_bytes(b"broken cache")
            corrupt = RemoteTabix(URL, directory, offline=True)
            with self.assertRaisesRegex(RemoteTabixError, "checksum mismatch"):
                _ = corrupt.references

    def test_offline_missing_window_does_not_download(self):
        blobs, _, _, _ = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            list(source.fetch("chr1", 1, 1))
            offline = RemoteTabix(URL, directory, offline=True)
            with self.assertRaisesRegex(RemoteTabixError, "Offline cache miss"):
                list(offline.fetch("chr1", 40000, 40000))

    def test_changed_remote_version_and_cached_headers_rejected(self):
        blobs, _, _, _ = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            list(source.fetch("chr1", 1, 1))
            source.etag = '"fixture-v2"'
            with self.assertRaisesRegex(RemoteTabixError, "etag changed"):
                list(source.fetch("chr1", 40000, 40000))
            entry = source.manifest["ranges"][0]
            header_path = source.cache_dir / entry["headers"]
            header_path.write_bytes(header_path.read_bytes() + b"tampered")
            offline = RemoteTabix(URL, directory, offline=True)
            with self.assertRaisesRegex(RemoteTabixError, "header checksum mismatch"):
                _ = offline.references

    def test_retry_and_failed_response_is_not_cached(self):
        blobs, _, _, _ = fixture()
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            original = source._request
            failures = []

            def once(url, start, end):
                if not failures:
                    failures.append(True)
                    return b"HTTP/2 200\r\n\r\n", b"ignored range"
                return original(url, start, end)

            source._request = once
            source.retries = 2
            from unittest.mock import patch
            with patch("lib_remote_tabix.time.sleep"):
                self.assertEqual(source.references, ("chr1", "chr2"))
            self.assertEqual(len(source.manifest["ranges"]), 1)

    def test_small_parallel_ranges_resume_a_failed_batch(self):
        blobs, _, _, _ = fixture()
        blobs[URL] = bytes(range(256)) * 6144
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            self.assertEqual(source.max_span_bytes, 256 * 1024)
            original = source._request
            active = maximum = 0
            failed = []
            lock = threading.Lock()

            def slow(url, start, end):
                nonlocal active, maximum
                with lock:
                    active += 1
                    maximum = max(maximum, active)
                try:
                    time.sleep(0.02)
                    if start == 2 * source.max_span_bytes and not failed:
                        failed.append(True)
                        raise RemoteTabixError("simulated stalled range")
                    return original(url, start, end)
                finally:
                    with lock:
                        active -= 1

            source._request = slow
            with self.assertRaisesRegex(RemoteTabixError, "simulated stalled range"):
                source._read_range(URL, 0, len(blobs[URL]) - 1)
            self.assertEqual(maximum, 4)
            saved = {(r['start'], r['end']) for r in source.manifest['ranges']}
            self.assertGreater(len(saved), 1)
            self.assertEqual(source._read_range(URL, 0, len(blobs[URL]) - 1), blobs[URL])
            # Successful parallel siblings were not requested again after failure.
            for interval in saved:
                self.assertEqual(sum((start, end) == interval for _, start, end in source.requests), 1)
            self.assertTrue(all(end - start + 1 <= 256 * 1024 for _, start, end in source.requests))
            offline = RemoteTabix(URL, directory, offline=True)
            self.assertEqual(offline._read_range(URL, 0, len(blobs[URL]) - 1), blobs[URL])

    def test_process_wide_network_limit_and_stall_cutoff(self):
        from unittest.mock import patch
        active = maximum = 0
        commands = []
        lock = threading.Lock()

        def network(command, **kwargs):
            nonlocal active, maximum
            with lock:
                active += 1
                maximum = max(maximum, active)
                commands.append(command)
            try:
                time.sleep(0.04)
                return SimpleNamespace(returncode=7, stderr=b"simulated network failure")
            finally:
                with lock:
                    active -= 1

        with tempfile.TemporaryDirectory() as directory:
            sources = [RemoteTabix(URL + str(i), directory) for i in range(8)]

            def call(source):
                with self.assertRaises(RemoteTabixError):
                    source._request(source.url, 0, 1)

            with patch("lib_remote_tabix.subprocess.run", side_effect=network):
                with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                    list(pool.map(call, sources))
            self.assertEqual(maximum, 4)
            self.assertEqual(len(commands), 8)
            for command in commands:
                self.assertEqual(command[command.index('--max-time') + 1], '60')
                self.assertEqual(command[command.index('--speed-limit') + 1], '1024')
                self.assertEqual(command[command.index('--speed-time') + 1], '20')

    def test_range_validation_rejects_full_body_and_wrong_coordinates(self):
        samples = [
            (b"HTTP/1.1 200 OK\r\nContent-Length: 3\r\n\r\n", b"abc"),
            (b"HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 1-3/10\r\n\r\n", b"abc"),
            (b"HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-2/10\r\n\r\n", b"ab"),
            (b"HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-2/10\r\nContent-Encoding: gzip\r\n\r\n", b"abc"),
        ]
        for headers, body in samples:
            with self.assertRaises(RemoteTabixError):
                validate_range_response(headers, body, 0, 2)
        self.assertEqual(validate_range_response(
            b"HTTP/2 206\r\nContent-Range: bytes 0-2/3\r\n\r\n", b"abc", 0, 65535), (2, 3))

    def test_content_length_validation(self):
        header = b"HTTP/2 206\r\nContent-Range: bytes 0-2/3\r\nContent-Length: "
        invalid = (b"", b"abc", b"+3", b"-3", b"0_3", b"3.0", b"3, 3",
                   b"\xb3", b"0", b"2", b"4", b"9" * 5000)
        for length in invalid:
            with self.subTest(length=length[:20]):
                with self.assertRaisesRegex(RemoteTabixError, "Content-Length"):
                    validate_range_response(header + length + b"\r\n\r\n", b"abc", 0, 2)
        for length in (b"3", b"0003", b"0" * 5000 + b"3"):
            with self.subTest(valid_length=length[:20]):
                self.assertEqual(validate_range_response(
                    header + length + b"\r\n\r\n", b"abc", 0, 2), (2, 3))

    def test_invalid_content_length_retries_without_caching(self):
        from unittest.mock import patch
        header = b"HTTP/2 206\r\nContent-Range: bytes 0-2/3\r\nContent-Length: "
        valid = header + b"3\r\n\r\n"
        for length in (b"bad", b"2"):
            with self.subTest(length=length), tempfile.TemporaryDirectory() as directory:
                source = RemoteTabix(URL, directory, retries=2)
                invalid = header + length + b"\r\n\r\n"
                with patch.object(source, "_request", return_value=(invalid, b"abc")) as request, \
                        patch("lib_remote_tabix.time.sleep") as sleep:
                    with self.assertRaisesRegex(RemoteTabixError, "Failed after 2 range attempts"):
                        source._download(URL, 0, 2)
                    self.assertEqual(request.call_count, 2)
                    sleep.assert_called_once_with(1)
                self.assertEqual(source.manifest["ranges"], [])
                self.assertEqual(source.manifest["totals"], {})
                self.assertEqual(list(source.cache_dir.iterdir()), [])
                with patch.object(source, "_request", side_effect=[
                        (invalid, b"abc"), (valid, b"abc")]) as request, \
                        patch("lib_remote_tabix.time.sleep") as sleep:
                    entry = source._download(URL, 0, 2)
                    self.assertEqual(request.call_count, 2)
                    sleep.assert_called_once_with(1)
                self.assertEqual(len(source.manifest["ranges"]), 1)
                self.assertEqual(source._cached_bytes(entry), b"abc")
                self.assertEqual((source.cache_dir / entry["headers"]).read_bytes(), valid)

    def test_bgzf_crc_member_boundaries_and_virtual_ends(self):
        one, two = bgzf(b"first\n"), bgzf(b"second\n")
        self.assertEqual(list(bgzf_blocks(one + two)), [(0, b"first\n"), (len(one), b"second\n")])
        damaged = bytearray(one)
        damaged[-8] ^= 1
        with self.assertRaisesRegex(RemoteTabixError, "CRC"):
            list(bgzf_blocks(damaged))
        with self.assertRaisesRegex(RemoteTabixError, "Truncated"):
            list(bgzf_blocks(one[:-1]))
        joined = bytearray(one + two)
        struct.pack_into("<H", joined, 16, len(joined) - 1)
        with self.assertRaisesRegex(RemoteTabixError, "boundary"):
            list(bgzf_blocks(joined))
        blobs, _, _, _ = fixture()
        blobs[URL] = one + two + bgzf(b"")
        with tempfile.TemporaryDirectory() as directory:
            source = MemoryRemote(directory, blobs)
            self.assertEqual(source._chunk_data(0, len(one) << 16), b"first\n")
            self.assertEqual(source._chunk_data(2, (len(one) << 16) | 3), b"rst\nsec")
            with self.assertRaisesRegex(RemoteTabixError, "Virtual offset"):
                source._chunk_data(0, (len(one) << 16) | 20)

    def test_uint64_offsets_and_bad_indexes(self):
        offset = (1 << 56) + 3
        names = b"chr1\0"
        raw = (b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 2, 2, 35, 0, len(names)) + names
               + struct.pack("<iIiQQiQ", 1, 4681, 1, offset, offset + 100, 1, offset))
        index = TabixIndex(raw)
        self.assertEqual(index.chunks("chr1", 0, 1), [(offset, offset + 100)])
        with self.assertRaises(RemoteTabixError):
            TabixIndex(raw[:-1])
        with self.assertRaises(RemoteTabixError):
            TabixIndex(b"CSI\x01" + raw[4:])
        with self.assertRaises(ValueError):
            reg2bins(0, (1 << 29) + 1)


if __name__ == "__main__":
    unittest.main()
