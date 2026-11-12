"""Offline checks for safe selective TAR access and phenotype SD extraction."""

import gzip
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/fetch_gtex_expression.py"
TEST_DIRECTORY = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("expression_fetch", MODULE_PATH)
fetch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fetch)


class MemoryRanges:
    def __init__(self, archive, directory):
        self.archive = archive
        self.directory = directory
        self.total_size = len(archive)
        self.ranges = []

    def get(self, start, length, purpose):
        self.ranges.append((start, length, purpose))
        path = self.directory / f"{start}-{length}.bin"
        path.write_bytes(self.archive[start:start + length])
        return path


class ExpressionAcquisitionTests(unittest.TestCase):
    def make_archive(self, format):
        raw = io.BytesIO()
        with tarfile.open(fileobj=raw, mode="w", format=format) as archive:
            # Long path exercises both PAX and GNU long-name extension handling.
            prefix = "long_directory_" * 12
            for tissue in fetch.TISSUES:
                bed = (
                    "#chr\tstart\tend\tgene_id\tsample1\tsample2\tsample3\n"
                    "chr1\t99\t100\tENSG00000000001.1\t1\t2\t3\n"
                    "chr1\t199\t200\tENSG00000000002.2\t0\t1\t4\n"
                ).encode()
                data = gzip.compress(bed, mtime=0)
                member = tarfile.TarInfo(f"{prefix}/{tissue}.v8.normalized_expression.bed.gz")
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
            member = tarfile.TarInfo("large_unrelated_member")
            member.size = 100_000
            archive.addfile(member, io.BytesIO(b"x" * member.size))
        return raw.getvalue()

    def test_pax_and_gnu_access_and_sd(self):
        for archive_format in (tarfile.PAX_FORMAT, tarfile.GNU_FORMAT):
            with self.subTest(format=archive_format), tempfile.TemporaryDirectory(dir=TEST_DIRECTORY) as temporary:
                directory = Path(temporary)
                self.assertTrue(directory.resolve().is_relative_to(TEST_DIRECTORY))
                archive = self.make_archive(archive_format)
                cache = MemoryRanges(archive, directory)
                members = fetch.fetch_members(cache, directory)
                rows, coverage = fetch.extract_sd(members, {"ENSG00000000001": "ENSG00000000001.1"})
                self.assertEqual(len(rows), 2)
                self.assertTrue(all(row["n"] == 3 and row["sdY"] == "1" for row in rows))
                self.assertTrue(all(row["tss"] == 100 for row in rows))
                self.assertTrue(all(row["observed_count"] == 1 for row in coverage))
                self.assertLess(sum(length for _, length, _ in cache.ranges), len(archive) / 10)
                self.assertEqual(cache.ranges[-1][2], "expression:Colon_Transverse")

    def test_bad_tar_checksum_is_rejected_before_payload(self):
        with tempfile.TemporaryDirectory(dir=TEST_DIRECTORY) as temporary:
            self.assertTrue(Path(temporary).resolve().is_relative_to(TEST_DIRECTORY))
            archive = bytearray(self.make_archive(tarfile.PAX_FORMAT))
            archive[0] ^= 1
            cache = MemoryRanges(archive, Path(temporary))
            with self.assertRaises(tarfile.InvalidHeaderError):
                fetch.fetch_members(cache, Path(temporary))
            self.assertEqual(len(cache.ranges), 1)

    def test_http_range_integrity(self):
        headers = "HTTP/1.1 200 Connection established\r\n\r\nHTTP/2 206\r\nContent-Range: bytes 512-1023/2048\r\nContent-Length: 512\r\n\r\n"
        total, _ = fetch.checked_content_range(headers, 512, 1023)
        self.assertEqual(total, 2048)
        with self.assertRaises(ValueError):
            fetch.checked_content_range(headers.replace("HTTP/2 206", "HTTP/2 200"), 512, 1023)
        with self.assertRaises(ValueError):
            fetch.checked_content_range(headers, 0, 511)

    def test_content_length_uses_ascii_decimal_notation(self):
        header = "HTTP/2 206\r\nContent-Range: bytes 0-2/3\r\n"
        invalid = ("", "abc", "+3", "-3", "0_3", "3.0", "3, 3",
                   "\u0663", "\xb3", "0", "2", "4", "9" * 5000)
        for length in invalid:
            with self.subTest(length=length[:20]):
                with self.assertRaisesRegex(ValueError, "Content-Length"):
                    fetch.checked_content_range(
                        header + f"Content-Length: {length}\r\n\r\n", 0, 2)
        for length in ("3", "0003", "0" * 5000 + "3"):
            with self.subTest(valid_length=length[:20]):
                total, _ = fetch.checked_content_range(
                    header + f"Content-Length: {length}\r\n\r\n", 0, 2)
                self.assertEqual(total, 3)
        self.assertEqual(fetch.checked_content_range(header + "\r\n", 0, 2)[0], 3)

    def test_range_cache_validates_length_before_commit_and_offline_reuse(self):
        from unittest.mock import patch
        url = "https://example.invalid/expression.tar"
        with tempfile.TemporaryDirectory(dir=TEST_DIRECTORY) as temporary:
            directory = Path(temporary)
            cache = fetch.RangeCache(url, directory, offline=False, curl="curl.exe")
            response_length = "+3"

            def download(command, **kwargs):
                Path(command[command.index("--output") + 1]).write_bytes(b"abc")
                headers = ("HTTP/2 206\r\nContent-Range: bytes 0-2/3\r\n"
                           f"Content-Length: {response_length}\r\n\r\n")
                Path(command[command.index("--dump-header") + 1]).write_bytes(headers.encode("ascii"))

            with patch.object(fetch.subprocess, "run", side_effect=download) as request:
                with self.assertRaisesRegex(ValueError, "Content-Length"):
                    cache.get(0, 3, "fixture")
                self.assertEqual(cache.used, [])
                self.assertIsNone(cache.total_size)
                self.assertIsNone(cache.etag)
                for pattern in ("*.bin", "*.headers.txt", "*.json"):
                    self.assertEqual(list(directory.glob(pattern)), [])
                response_length = "0" * 5000 + "3"
                self.assertEqual(cache.get(0, 3, "fixture").read_bytes(), b"abc")
                self.assertEqual(request.call_count, 2)
                self.assertEqual(len(cache.used), 1)
            offline = fetch.RangeCache(url, directory, offline=True, curl="curl.exe")
            with patch.object(fetch.subprocess, "run") as request:
                self.assertEqual(offline.get(0, 3, "fixture").read_bytes(), b"abc")
                request.assert_not_called()
            self.assertEqual(len(offline.used), 1)

    def test_gene_version_mismatch_is_not_silently_mapped(self):
        with tempfile.TemporaryDirectory(dir=TEST_DIRECTORY) as temporary:
            directory = Path(temporary)
            self.assertTrue(directory.resolve().is_relative_to(TEST_DIRECTORY))
            cache = MemoryRanges(self.make_archive(tarfile.PAX_FORMAT), directory)
            members = fetch.fetch_members(cache, directory)
            with self.assertRaisesRegex(ValueError, "version mismatch"):
                fetch.extract_sd(members, {"ENSG00000000001": "ENSG00000000001.2"})


if __name__ == "__main__":
    unittest.main()
