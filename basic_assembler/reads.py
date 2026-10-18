"""Validate FASTQ input and build one reusable, strand-independent read pool.

Subsampling takes the first N physical reads across all single-end files, or
the first N pairs across all paired lanes. Only consumed records are validated.
Optional prepared FASTQ contains exactly the retained, trimmed physical reads;
reverse complements are added to the in-memory pool, never to that file.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, field
import gzip
import math
import os
from pathlib import Path
import re
from typing import Iterator, TextIO


_IUPAC = frozenset("ACGTRYSWKMBDHVN")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_AMBIGUOUS = re.compile("[RYSWKMBDHV]")
_CASAVA_MATE = re.compile(r"^[12]:[YN]:[0-9]+:")
_INVALID_QUALITY = re.compile(r"[^!-~]")


class ReadValidationError(ValueError):
    """An input file or record cannot be interpreted safely."""


@dataclass
class ReadPool:
    counts: dict[str, int]
    quality_weights: dict[str, float]
    mate_links: dict[str, dict[str, int]]
    max_read_length: int
    stats: dict
    ambiguous_sequences: dict[str, str] = field(default_factory=dict)


def _revcom(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@contextmanager
def _open_fastq(path: str | Path) -> Iterator[TextIO]:
    """Recognize gzip by its magic bytes, including files without a .gz suffix."""
    path = Path(path)
    try:
        with path.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        opener = gzip.open if compressed or path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="ascii", newline=None) as handle:
            yield handle
    except (OSError, EOFError, UnicodeError) as exc:
        raise ReadValidationError(f"Cannot read FASTQ {path}: {exc}") from exc


def _read_record(handle: TextIO, path: str | Path, number: int):
    header = handle.readline()
    if not header:
        return None
    lines = [header, handle.readline(), handle.readline(), handle.readline()]
    location = f"{path}, record {number}"
    if any(line == "" for line in lines[1:]):
        raise ReadValidationError(f"Truncated FASTQ at {location}: expected four lines")
    header, sequence, plus, quality = [line.rstrip("\r\n") for line in lines]
    if not header.startswith("@") or not header[1:].strip():
        raise ReadValidationError(f"Invalid FASTQ header at {location}")
    if not plus.startswith("+"):
        raise ReadValidationError(f"Invalid FASTQ separator at {location}")
    header_parts = header[1:].split()
    read_id = header_parts[0]
    mate = int(read_id[-1]) if read_id.endswith(("/1", "/2")) else None
    if len(header_parts) > 1 and _CASAVA_MATE.match(header_parts[1]):
        casava_mate = int(header_parts[1][0])
        if mate is not None and mate != casava_mate:
            raise ReadValidationError(f"Conflicting mate labels at {location}")
        mate = casava_mate
    if plus[1:].strip() and plus[1:].split()[0] != read_id:
        raise ReadValidationError(f"FASTQ separator identifier differs at {location}")
    sequence = sequence.upper().replace("U", "T")
    if not sequence or not set(sequence) <= _IUPAC:
        raise ReadValidationError(f"Empty or non-IUPAC sequence at {location}")
    if len(sequence) != len(quality):
        raise ReadValidationError(f"Sequence/quality length mismatch at {location}")
    if _INVALID_QUALITY.search(quality):
        raise ReadValidationError(f"Invalid Phred+33 quality at {location}")
    return read_id, sequence, quality, mate


def _mate_id(identifier: str) -> str:
    return identifier[:-2] if identifier.endswith(("/1", "/2")) else identifier


def _same_file(first: str | Path, second: str | Path) -> bool:
    if Path(first).resolve() == Path(second).resolve():
        return True
    try:
        return os.path.samefile(first, second)
    except OSError:
        return False


def load_read_pool(
    single: list[str],
    left: list[str],
    right: list[str],
    subsample: int | None = None,
    min_read_length: int = 16,
    trim_quality: int = 0,
    quality_weighted: bool = False,
    prepared_path: str | Path | None = None,
    retain_mates: bool = False,
) -> ReadPool:
    """Read FASTQ once, validating synchronized mates and retaining useful reads.

    Terminal Ns and low-quality terminal bases are trimmed together on both ends.
    Reads shorter than ``min_read_length`` or with >=50% one A/C/G/T nucleotide
    are omitted. IUPAC ambiguity is retained; U is normalized to T. Weighting,
    when enabled, sums the probability that every retained base is correct.
    A palindromic read contributes once, rather than twice, to its pool entry.

    With ``retain_mates``, all orientation combinations of a physical pair are
    linked; each directed edge receives at most one observation per pair. The
    default avoids storing these links. No original read IDs are retained.
    Prepared reads receive monotonically increasing synthetic IDs.
    Only prepared reads containing non-ACGTN IUPAC symbols are retained in the
    sparse ``ambiguous_sequences`` map so SAM output can restore their symbols.
    """
    single, left, right = list(single), list(left), list(right)
    if single and (left or right):
        raise ReadValidationError("Choose single-end input or paired-end input, not both")
    if not single and not (left and right):
        raise ReadValidationError("Provide single-end files or both paired-end file lists")
    if len(left) != len(right):
        raise ReadValidationError("Paired input must contain the same number of left and right lanes")
    if subsample is not None and (isinstance(subsample, bool) or not isinstance(subsample, int) or subsample < 1):
        raise ReadValidationError("subsample must be a positive integer or None")
    if isinstance(min_read_length, bool) or not isinstance(min_read_length, int) or min_read_length < 1:
        raise ReadValidationError("min_read_length must be a positive integer")
    if isinstance(trim_quality, bool) or not isinstance(trim_quality, int) or not 0 <= trim_quality <= 93:
        raise ReadValidationError("trim_quality must be an integer from 0 to 93")
    for first, second in zip(left, right):
        if _same_file(first, second):
            raise ReadValidationError("Paired inputs must not refer to the same file")
    if prepared_path is not None:
        destination = Path(prepared_path).resolve()
        if any(_same_file(source, destination) for source in single + left + right):
            raise ReadValidationError("Prepared FASTQ must not overwrite an input file")

    counts: Counter = Counter()
    weights: defaultdict = defaultdict(float)
    links: defaultdict = defaultdict(Counter)
    ambiguous_sequences = {}
    stats = {
        "input_mode": "single" if single else "paired",
        "input_files": len(single) + len(left) + len(right),
        "subsample_limit": subsample,
        "subsample_unit": "reads" if single else "pairs",
        "consumed_reads": 0,
        "consumed_pairs": 0,
        "retained_reads": 0,
        "retained_pairs": 0,
        "trimmed_reads": 0,
        "filtered_short_reads": 0,
        "filtered_low_complexity_reads": 0,
        "input_bases": 0,
        "retained_bases": 0,
        "mate_links_retained": retain_mates,
        "quality_weighted": quality_weighted,
        "ambiguous_read_sequences": 0,
    }
    maximum = 0

    def retain(record, prepared):
        nonlocal maximum
        _, sequence, quality, _ = record
        stats["consumed_reads"] += 1
        stats["input_bases"] += len(sequence)
        start, end = 0, len(sequence)
        while start < end and (sequence[start] == "N" or ord(quality[start]) - 33 < trim_quality):
            start += 1
        while end > start and (sequence[end - 1] == "N" or ord(quality[end - 1]) - 33 < trim_quality):
            end -= 1
        if start or end != len(sequence):
            stats["trimmed_reads"] += 1
        sequence, quality = sequence[start:end], quality[start:end]
        if len(sequence) < min_read_length:
            stats["filtered_short_reads"] += 1
            return None
        if any(sequence.count(base) * 2 >= len(sequence) for base in "ACGT"):
            stats["filtered_low_complexity_reads"] += 1
            return None
        stats["retained_reads"] += 1
        stats["retained_bases"] += len(sequence)
        maximum = max(maximum, len(sequence))
        if _AMBIGUOUS.search(sequence):
            stats["ambiguous_read_sequences"] += 1
            if prepared is not None:
                ambiguous_sequences[f"basic_{stats['consumed_reads']}"] = sequence
        oriented = {sequence, _revcom(sequence)}
        probability = math.prod(1.0 - 10.0 ** (-(ord(char) - 33) / 10.0)
                                for char in quality) if quality_weighted else 0.0
        for value in oriented:
            counts[value] += 1
            if quality_weighted:
                weights[value] += probability
        if prepared is not None:
            prepared.write(f"@basic_{stats['consumed_reads']}\n{sequence}\n+\n{quality}\n")
        return oriented

    with ExitStack() as stack:
        prepared = None
        if prepared_path is not None:
            prepared = stack.enter_context(open(prepared_path, "w", encoding="ascii", newline="\n"))
        if single:
            for source in single:
                if subsample is not None and stats["consumed_reads"] >= subsample:
                    break
                with _open_fastq(source) as handle:
                    number = 0
                    while subsample is None or stats["consumed_reads"] < subsample:
                        number += 1
                        record = _read_record(handle, source, number)
                        if record is None:
                            break
                        retain(record, prepared)
        else:
            for left_path, right_path in zip(left, right):
                if subsample is not None and stats["consumed_pairs"] >= subsample:
                    break
                with _open_fastq(left_path) as lhandle, _open_fastq(right_path) as rhandle:
                    number = 0
                    while subsample is None or stats["consumed_pairs"] < subsample:
                        number += 1
                        first = _read_record(lhandle, left_path, number)
                        second = _read_record(rhandle, right_path, number)
                        if first is None and second is None:
                            break
                        if first is None or second is None:
                            raise ReadValidationError(
                                f"Mate count mismatch in {left_path} and {right_path}, pair {number}")
                        if _mate_id(first[0]) != _mate_id(second[0]):
                            raise ReadValidationError(
                                f"Mate identifier mismatch at pair {number}: {first[0]!r} != {second[0]!r}")
                        if first[3] not in (None, 1) or second[3] not in (None, 2):
                            raise ReadValidationError(
                                f"Mate direction mismatch at pair {number}: left must be mate 1 and right mate 2")
                        stats["consumed_pairs"] += 1
                        a, b = retain(first, prepared), retain(second, prepared)
                        if a is not None and b is not None:
                            stats["retained_pairs"] += 1
                            if retain_mates:
                                edges = {(x, y) for x in a for y in b}
                                edges.update((y, x) for x in a for y in b)
                                for x, y in edges:
                                    links[x][y] += 1
    stats["unique_oriented_sequences"] = len(counts)
    return ReadPool(dict(counts), dict(weights),
                    {key: dict(value) for key, value in links.items()}, maximum, stats,
                    ambiguous_sequences)
