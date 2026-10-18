"""Exercise real subprocess handling with tiny executable alignment stubs."""

from pathlib import Path
import os
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import pytest

from basic_assembler.alignment import (AlignmentError, resolve_bowtie,
                                      run_alignment, validate_indexes)


SEQ = "ACGTACGTTGCACTGA"


def sam(identifier="read", reference="X|IGHV1*01|human", sequence=SEQ,
        mapq=42, flag=0, cigar="16M"):
    return f"{identifier}\t{flag}\t{reference}\t1\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t{'I' * len(sequence)}\n"


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="basic-align-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def executable(self, stdout="", stderr="", code=0):
        directory = self.root / "binary with spaces"
        directory.mkdir(exist_ok=True)
        path = directory / "bowtie2"
        path.write_text(f"#!{sys.executable}\nimport sys\n"
                        "if '--version' in sys.argv:\n print('bowtie2 version test'); sys.exit(0)\n"
                        f"sys.stdout.write({stdout!r})\nsys.stderr.write({stderr!r})\nsys.exit({code})\n")
        path.chmod(0o755)
        return str(path)

    def align(self, executable, **kwargs):
        return run_alignment(executable, self.root / "db with spaces" / "hv",
                             [str(self.root / "reads with spaces.fq")], 2,
                             self.root / "logs" / "align.stderr", **kwargs)

    def test_spaces_are_preserved_and_version_checked(self):
        executable = self.executable(sam())
        resolved, version = resolve_bowtie(Path(executable).parent)
        self.assertEqual(resolved, str(Path(executable).resolve()))
        self.assertIn("version test", version)
        result = self.align(executable, strand="v")
        self.assertEqual(result.mapped_reads, 1)
        self.assertIn("--norc", result.command)
        self.assertEqual(result.command[-1], str(self.root / "reads with spaces.fq"))
        self.assertEqual(result.candidates[0]["locus"], "IGH")

    def test_failure_is_propagated_even_after_valid_output(self):
        executable = self.executable(sam(), "index is broken\n", 7)
        with self.assertRaisesRegex(AlignmentError, "exit 7") as caught:
            self.align(executable)
        self.assertIn("index is broken", str(caught.exception))
        self.assertIn("index is broken", (self.root / "logs/align.stderr").read_text())

    def test_malformed_sam_is_rejected(self):
        values = ["incomplete\tSAM\n", sam(mapq="bad"), sam(mapq=300),
                  sam(sequence="AC?T"), sam(cigar="15M"), sam(cigar="*"),
                  sam(reference="*"), sam(reference=""), "@not-a-SAM-header\n"]
        for value in values:
            with self.subTest(value=value), self.assertRaises(AlignmentError):
                self.align(self.executable(value))

    def test_grouping_and_ranking_are_deterministic(self):
        rows = [sam("a", "Z|TRBV1*01|human", mapq=20),
                sam("b", "A|TRAV1*01|human", mapq=10),
                sam("c", "Z|TRBV1*01|human", mapq=40),
                sam("secondary", "Z|TRBV1*01|human", flag=256),
                sam("supplementary", "Z|TRBV1*01|human", flag=2048)]
        first = self.align(self.executable("".join(rows)))
        second = self.align(self.executable("".join(reversed(rows))))
        self.assertEqual(first.candidates, second.candidates)
        self.assertEqual(first.mapped_reads, 3)
        self.assertEqual(first.candidates[0]["count"], 2)
        self.assertEqual(first.candidates[0]["mapq"], 40)
        self.assertEqual(first.candidates[0]["locus"], "TRB")

    def test_empty_alignment_and_unmapped_records(self):
        self.assertEqual(self.align(self.executable()).candidates, [])
        unmapped = "x\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
        self.assertEqual(self.align(self.executable(unmapped)).mapped_reads, 0)

    def test_local_seeding_preserves_full_soft_clipped_sequence(self):
        result = self.align(self.executable(sam(cigar="5S11M")), local=True)
        self.assertIn("--very-sensitive-local", result.command)
        self.assertNotIn("--very-sensitive", result.command)
        self.assertEqual(result.command[result.command.index("-L") + 1], "10")
        self.assertEqual(result.command[result.command.index("-N") + 1], "1")
        self.assertEqual(result.candidates[0]["sequence"], SEQ)

    def test_iupac_restoration_respects_sam_strand_and_validates_alignment(self):
        original = SEQ[:5] + "Y" + SEQ[6:]
        emitted = original.replace("Y", "N")
        complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
        reverse = original.translate(complement)[::-1]
        rows = sam("basic_1", sequence=emitted) + sam("basic_2", sequence=emitted.translate(complement)[::-1], flag=16)
        result = self.align(self.executable(rows), original_sequences={"basic_1": original, "basic_2": original})
        self.assertEqual({item["sequence"] for item in result.candidates}, {original, reverse})
        self.assertEqual(result.mapped_reads, 2)
        wrong = "T" + emitted[1:]
        with self.assertRaisesRegex(AlignmentError, "SAM/original sequence mismatch"):
            self.align(self.executable(sam("basic_1", sequence=wrong)), original_sequences={"basic_1": original})

    def test_small_large_and_incomplete_indexes(self):
        for region in ("hv", "hc", "lv", "lc"):
            extension = "bt2l" if region == "lv" else "bt2"
            for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
                (self.root / f"{region}.{suffix}.{extension}").write_bytes(b"index")
        self.assertEqual(set(validate_indexes(self.root)), {"hv", "hc", "lv", "lc"})
        (self.root / "hc.rev.2.bt2").unlink()
        with self.assertRaisesRegex(AlignmentError, "hc"):
            validate_indexes(self.root)

    def test_nonexistent_executable(self):
        with self.assertRaises(AlignmentError):
            resolve_bowtie(self.root / "missing-bowtie")

    def assert_child_stopped(self, marker):
        # The child writes this marker after a delay unless it was stopped.
        # This needs no global process inspection.
        time.sleep(1.1)
        self.assertFalse(marker.exists(), "Aligner descendant survived cleanup")

    def child_spawning_executable(self, body):
        path = Path(self.executable())
        marker = self.root / "child-survived"
        ready = self.root / "child-ready"
        child = ("import time; from pathlib import Path; "
                 f"Path({str(ready)!r}).write_text('ready'); time.sleep(1.0); "
                 f"Path({str(marker)!r}).write_text('survived')")
        path.write_text(f"#!{sys.executable}\nimport subprocess, sys, time\nfrom pathlib import Path\n"
                        f"child = subprocess.Popen([sys.executable, '-c', {child!r}])\n"
                        f"while not Path({str(ready)!r}).exists(): time.sleep(0.005)\n" + body)
        return str(path), marker

    @unittest.skipUnless(os.name == "posix", "Process-group cleanup is POSIX-specific")
    def test_timeout_stops_descendant_holding_stdout_after_wrapper_exits(self):
        executable, marker = self.child_spawning_executable(
            "print('waiting for child', file=sys.stderr, flush=True)\n")
        started = time.monotonic()
        with self.assertRaisesRegex(AlignmentError, "timed out") as error:
            self.align(executable, timeout=0.5)
        self.assertLess(time.monotonic() - started, 3)
        # Include whatever the child emitted; startup may itself exhaust the
        # timeout on heavily instrumented CI workers.
        diagnostics = (self.root / "logs/align.stderr").read_text().strip()
        self.assertIn(diagnostics, str(error.exception))
        self.assert_child_stopped(marker)

    @unittest.skipUnless(os.name == "posix", "Process-group cleanup is POSIX-specific")
    def test_bad_sam_stops_descendants(self):
        executable, marker = self.child_spawning_executable(
            "print('bad SAM', flush=True)\ntime.sleep(60)\n")
        with self.assertRaisesRegex(AlignmentError, "Malformed SAM"):
            self.align(executable, timeout=5)
        self.assert_child_stopped(marker)

    @unittest.skipUnless(os.name == "posix", "Process-group cleanup is POSIX-specific")
    def test_version_timeout_stops_descendants(self):
        executable, marker = self.child_spawning_executable("time.sleep(60)\n")
        with patch("basic_assembler.alignment._VERSION_TIMEOUT", 0.5):
            with self.assertRaisesRegex(AlignmentError, "version query timed out"):
                resolve_bowtie(executable)
        self.assert_child_stopped(marker)

    def test_invalid_timeouts_and_successful_bounded_alignment(self):
        executable = self.executable(sam())
        for timeout in (0, -1, float("inf"), float("nan"), True, "bad"):
            with self.subTest(timeout=timeout), self.assertRaises(ValueError):
                self.align(executable, timeout=timeout)
        self.assertEqual(self.align(executable, timeout=10).mapped_reads, 1)


if __name__ == "__main__":
    unittest.main()


@pytest.mark.integration
def test_real_bowtie_iupac_symbols_are_preserved_in_both_input_orientations(tmp_path, repository_root, real_bowtie2):
    from basic_assembler.assembly import revcom
    from basic_assembler.reads import load_read_pool
    from benchmarks.generate_synthetic import _references
    sequence = _references(repository_root / "db/human/BCR/hv")[0][1][30:105]
    symbol = {"A": "R", "G": "R", "C": "Y", "T": "Y"}[sequence[35]]
    sequence = sequence[:35] + symbol + sequence[36:]
    source, prepared = tmp_path / "reads.fq", tmp_path / "prepared.fq"
    source.write_text("".join(f"@read{i}\n{value}\n+\n{'I' * len(value)}\n"
                              for i, value in enumerate((sequence, revcom(sequence)))))
    pool = load_read_pool([str(source)], [], [], prepared_path=prepared)
    result = run_alignment(real_bowtie2, repository_root / "db/human/BCR/hv",
                           [str(prepared)], 1, tmp_path / "alignment.stderr",
                           original_sequences=pool.ambiguous_sequences)
    assert result.mapped_reads == 2
    assert sum(item["count"] for item in result.candidates) == 2
    assert all(item["sequence"] == sequence and item["sequence"] in pool.counts
               for item in result.candidates)
