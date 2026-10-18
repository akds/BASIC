"""FASTQ preparation regression tests; also runnable with unittest."""

import gzip
import math
import os
from pathlib import Path
import tempfile
import unittest

from basic_assembler.reads import ReadValidationError, load_read_pool


SEQ = "ACGTACGTTGCACTGA"
MATE = "TGCAAGTCCGATGTAC"


def revcom(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class ReadPoolTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="basic-reads-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fastq(self, name, records, compressed=False):
        path = self.root / name
        text = "".join(f"@{identifier}\n{sequence}\n+\n{quality}\n"
                       for identifier, sequence, quality in records)
        if compressed:
            path.write_bytes(gzip.compress(text.encode("ascii")))
        else:
            path.write_text(text, encoding="ascii")
        return str(path)

    def test_orientation_palindrome_and_prepared_reads(self):
        palindrome = "ACGT" * 4
        path = self.fastq("reads.fq", [("x", SEQ.lower(), "I" * 16),
                                      ("y", palindrome, "I" * 16)])
        prepared = self.root / "prepared.fq"
        pool = load_read_pool([path], [], [], prepared_path=prepared)
        self.assertEqual(pool.counts, {SEQ: 1, revcom(SEQ): 1, palindrome: 1})
        self.assertEqual(pool.max_read_length, 16)
        self.assertEqual(pool.stats["retained_reads"], 2)
        self.assertEqual(prepared.read_text().count("@basic_"), 2)
        self.assertEqual(pool.quality_weights, {})

    def test_gzip_detected_by_magic_and_iupac_normalized(self):
        sequence = "ACGTRYSWKMBDHVAC"
        path = self.fastq("compressed.data", [("x", sequence.lower(), "I" * 16)], True)
        pool = load_read_pool([path], [], [])
        self.assertIn(sequence, pool.counts)

    def test_sparse_iupac_restoration_map_uses_prepared_identifiers(self):
        sequence = "ACGTRYSWKMBDHVAC"
        path = self.fastq("ambiguity.fq", [("original-1", SEQ, "I" * 16),
                                          ("original-2", sequence, "I" * 16)])
        prepared = self.root / "prepared.fq"
        pool = load_read_pool([path], [], [], prepared_path=prepared)
        self.assertEqual(pool.ambiguous_sequences, {"basic_2": sequence})
        self.assertEqual(pool.stats["ambiguous_read_sequences"], 1)
        self.assertIn(sequence, pool.counts)
        self.assertIn("@basic_2\n" + sequence, prepared.read_text())

    def test_trimming_filtering_and_quality_weights(self):
        path = self.fastq("trim.fq", [("x", "N" + SEQ + "NN", "!" + "I" * 16 + "!!"),
                                     ("short", "N" + SEQ[:-1], "I" * 16),
                                     ("low", "A" * 8 + "CGTACGTC", "I" * 16)])
        prepared = self.root / "out.fq"
        pool = load_read_pool([path], [], [], trim_quality=20,
                              quality_weighted=True, prepared_path=prepared)
        self.assertEqual(pool.stats["retained_reads"], 1)
        self.assertEqual(pool.stats["filtered_short_reads"], 1)
        self.assertEqual(pool.stats["filtered_low_complexity_reads"], 1)
        self.assertAlmostEqual(pool.quality_weights[SEQ], math.pow(0.9999, 16))
        self.assertEqual(prepared.read_text().splitlines()[1], SEQ)

    def test_terminal_ns_do_not_hide_low_quality_bases(self):
        path = self.fastq("hidden-low-quality.fq", [("x", "NA" + SEQ + "CN", "I!" + "I" * 16 + "!I")])
        prepared = self.root / "trimmed.fq"
        pool = load_read_pool([path], [], [], trim_quality=20, prepared_path=prepared)
        self.assertEqual(pool.counts, {SEQ: 1, revcom(SEQ): 1})
        self.assertEqual(prepared.read_text().splitlines()[1:4], [SEQ, "+", "I" * 16])

    def test_single_subsampling_counts_physical_reads_across_files(self):
        first = self.fastq("a.fq", [("x", "A" * 16, "I" * 16)])
        second = self.fastq("b.fq", [("y", SEQ, "I" * 16), ("z", MATE, "I" * 16)])
        pool = load_read_pool([first, second], [], [], subsample=2)
        self.assertEqual(pool.stats["consumed_reads"], 2)
        self.assertEqual(pool.stats["retained_reads"], 1)
        self.assertNotIn(MATE, pool.counts)

    def test_pair_subsampling_across_lanes_and_mate_orientations(self):
        left = [self.fastq(f"l{n}.fq", [(f"cell{n}/1", SEQ, "I" * 16)]) for n in range(3)]
        right = [self.fastq(f"r{n}.fq", [(f"cell{n}/2", MATE, "I" * 16)]) for n in range(3)]
        prepared = self.root / "pairs.fq"
        pool = load_read_pool([], left, right, subsample=2, prepared_path=prepared,
                              retain_mates=True)
        self.assertEqual(pool.stats["consumed_reads"], 4)
        self.assertEqual(pool.stats["consumed_pairs"], 2)
        self.assertEqual(pool.stats["retained_pairs"], 2)
        self.assertEqual(pool.mate_links[SEQ][revcom(MATE)], 2)
        self.assertEqual(pool.mate_links[MATE][SEQ], 2)
        self.assertEqual(prepared.read_text().count("@basic_"), 4)

    def test_pair_stats_without_retaining_links_and_filtered_orphans(self):
        left = self.fastq("left.fq", [("one/1", SEQ, "I" * 16),
                                     ("two/1", SEQ, "I" * 16)])
        right = self.fastq("right.fq", [("one/2", MATE, "I" * 16),
                                       ("two/2", "AAAA", "IIII")])
        pool = load_read_pool([], [left], [right])
        self.assertEqual(pool.mate_links, {})
        self.assertEqual(pool.quality_weights, {})
        self.assertEqual(pool.stats["consumed_pairs"], 2)
        self.assertEqual(pool.stats["retained_pairs"], 1)
        self.assertEqual(pool.stats["retained_reads"], 3)
        self.assertEqual(pool.stats["filtered_short_reads"], 1)

    def test_pair_id_and_count_mismatches(self):
        left = self.fastq("left.fq", [("one/1", SEQ, "I" * 16)])
        wrong = self.fastq("wrong.fq", [("two/2", MATE, "I" * 16)])
        empty = self.fastq("empty.fq", [])
        with self.assertRaisesRegex(ReadValidationError, "identifier mismatch"):
            load_read_pool([], [left], [wrong])
        with self.assertRaisesRegex(ReadValidationError, "count mismatch"):
            load_read_pool([], [left], [empty])
        with self.assertRaisesRegex(ReadValidationError, "same number"):
            load_read_pool([], [left, left], [wrong])

    def test_mate_directions_and_conflicting_labels_are_rejected(self):
        for first, second in [("x/2", "x/1"), ("x/1", "x/1"),
                              ("x 2:N:0:ACGT", "x 1:N:0:ACGT"),
                              ("x/1 2:N:0:ACGT", "x/2")]:
            with self.subTest(first=first, second=second):
                left = self.fastq("left.fq", [(first, SEQ, "I" * 16)])
                right = self.fastq("right.fq", [(second, MATE, "I" * 16)])
                with self.assertRaisesRegex(ReadValidationError, "[Mm]ate"):
                    load_read_pool([], [left], [right])
        # SRA ordinal metadata is not a CASAVA mate direction.
        left = self.fastq("left.fq", [("SRR1.1 1 length=16", SEQ, "I" * 16)])
        right = self.fastq("right.fq", [("SRR1.1 1 length=16", MATE, "I" * 16)])
        self.assertEqual(load_read_pool([], [left], [right]).stats["consumed_pairs"], 1)

    def test_malformed_fastq_is_rejected(self):
        bad_contents = ["@x\nACGT\n+\n", "x\nACGT\n+\nIIII\n",
                        "@x\nACGT\n-\nIIII\n", "@x\nACGT\n+\nIII\n",
                        "@x\nAC?T\n+\nIIII\n", "@x\nACGT\n+\nI I!\n",
                        "@x\nACGT\n+y\nIIII\n"]
        for content in bad_contents:
            with self.subTest(content=content):
                path = self.root / "bad.fq"
                path.write_text(content)
                with self.assertRaises(ReadValidationError):
                    load_read_pool([str(path)], [], [])

    def test_phred33_quality_range_includes_both_endpoints_only(self):
        for quality in ("!" * 16, "~" * 16):
            path = self.fastq("valid-quality.fq", [("x", SEQ, quality)])
            self.assertEqual(load_read_pool([path], [], []).stats["retained_reads"], 1)
        for character in (" ", "\x00", "\x7f"):
            path = self.fastq("bad-quality.fq", [("x", SEQ, "I" * 8 + character + "I" * 7)])
            with self.subTest(character=repr(character)), self.assertRaises(ReadValidationError):
                load_read_pool([path], [], [])

    def test_bad_gzip_and_overwriting_input_are_rejected(self):
        path = self.root / "bad.fq.gz"
        path.write_bytes(b"not gzip")
        with self.assertRaises(ReadValidationError):
            load_read_pool([str(path)], [], [])
        normal = self.fastq("normal.fq", [("x", SEQ, "I" * 16)])
        with self.assertRaisesRegex(ReadValidationError, "overwrite"):
            load_read_pool([normal], [], [], prepared_path=normal)

    def test_hardlinked_inputs_cannot_be_overwritten_or_used_as_two_mates(self):
        path = self.fastq("original.fq", [("x", SEQ, "I" * 16)])
        alias = self.root / "alias.fq"
        os.link(path, alias)
        original = alias.read_bytes()
        with self.assertRaisesRegex(ReadValidationError, "overwrite"):
            load_read_pool([path], [], [], prepared_path=alias)
        with self.assertRaisesRegex(ReadValidationError, "same file"):
            load_read_pool([], [path], [str(alias)])
        self.assertEqual(alias.read_bytes(), original)

    def test_invalid_options(self):
        path = self.fastq("reads.fq", [("x", SEQ, "I" * 16)])
        for kwargs in ({"subsample": 0}, {"subsample": -1}, {"subsample": True},
                       {"min_read_length": 0}, {"min_read_length": True},
                       {"trim_quality": 94}, {"trim_quality": False}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ReadValidationError):
                load_read_pool([path], [], [], **kwargs)
        with self.assertRaises(ReadValidationError):
            load_read_pool([path], [path], [path])


if __name__ == "__main__":
    unittest.main()
