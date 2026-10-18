"""Checked Bowtie2 execution and deterministic germline seed candidates."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import threading

from .processes import cleanup_process, stop_process_tree as _stop_process_tree


_VERSION_TIMEOUT = 30.0
_IUPAC = frozenset("ACGTRYSWKMBDHVN")
_AMBIGUOUS = frozenset("RYSWKMBDHV")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_INVALID_QUALITY = re.compile(r"[^!-~]")


class AlignmentError(RuntimeError):
    """Bowtie2 failed, or its output is not valid SAM."""


@dataclass
class MappingResult:
    candidates: list[dict]
    mapped_reads: int
    command: list[str]


def resolve_bowtie(executable: str | Path = "bowtie2") -> tuple[str, str]:
    """Resolve a directory, explicit executable, or executable in PATH."""
    requested = Path(executable).expanduser()
    if requested.is_dir():
        requested /= "bowtie2"
    resolved = str(requested.resolve()) if requested.is_file() else shutil.which(str(requested))
    if not resolved or not os.access(resolved, os.X_OK):
        raise AlignmentError(f"Bowtie2 executable not found or not executable: {executable}")
    try:
        process = subprocess.Popen([resolved, "--version"], stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True,
                                   start_new_session=(os.name == "posix"))
    except OSError as exc:
        raise AlignmentError(f"Cannot query Bowtie2 version at {resolved}: {exc}") from exc
    try:
        output, _ = process.communicate(timeout=_VERSION_TIMEOUT)
    except BaseException as exc:
        cleanup_process(process)
        if isinstance(exc, subprocess.TimeoutExpired):
            raise AlignmentError(f"Bowtie2 version query timed out after {_VERSION_TIMEOUT:g} seconds at {resolved}") from exc
        raise
    if process.returncode:
        raise AlignmentError(f"Bowtie2 --version failed with exit {process.returncode}: {output.strip()}")
    lines = output.strip().splitlines()
    if not lines:
        raise AlignmentError(f"Bowtie2 --version produced no version information: {resolved}")
    return resolved, lines[0]


def validate_indexes(database_path: str | Path) -> dict[str, str]:
    """Require all six readable, nonempty shards of each small or large index."""
    root = Path(database_path)
    suffixes = ("1", "2", "3", "4", "rev.1", "rev.2")
    indexes = {}
    for region in ("hv", "hc", "lv", "lc"):
        prefix = root / region
        valid = False
        for extension in ("bt2", "bt2l"):
            shards = [Path(f"{prefix}.{suffix}.{extension}") for suffix in suffixes]
            try:
                valid = all(shard.is_file() and os.access(shard, os.R_OK)
                            and shard.stat().st_size > 0 for shard in shards)
            except OSError:
                valid = False
            if valid:
                break
        if not valid:
            raise AlignmentError(
                f"Incomplete Bowtie2 index {prefix}: expected six readable, nonempty .bt2 or .bt2l shards")
        indexes[region] = str(prefix)
    return indexes


def _stderr_tail(path: Path) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            handle.seek(max(0, handle.tell() - 2000))
            return handle.read().decode("utf-8", errors="replace").strip()
    except OSError:
        return ""


def run_alignment(
    bowtie: str | Path,
    index_prefix: str | Path,
    read_paths: list[str],
    threads: int,
    stderr_path: str | Path,
    strand: str | None = None,
    local: bool = False,
    timeout: float | None = None,
    original_sequences: Mapping[str, str] | None = None,
) -> MappingResult:
    """Stream primary mapped SAM records and group by sequence and reference.

    Candidate ordering is descending read count, descending maximum MAPQ, then
    reference and sequence lexicographically. Counts describe physical mapped
    reads. SAM SEQ is already oriented to the reference; do not reverse it again.
    ``strand`` accepts v/--norc or c/--nofw (and the region names hv/hc/lv/lc).
    Local seeding uses a sensitive, short-seed preset for exploratory homology
    detection. It retains the complete SAM read; homology is not a validated
    species or gene assignment.
    When supplied, ``timeout`` bounds the entire alignment subprocess and its
    output stream. On POSIX, cleanup attempts to stop the process group, falling
    back to the direct child if group signalling is denied. No fixed limit is
    imposed on legitimate long-running jobs by default.
    ``original_sequences`` is an optional sparse prepared-QNAME lookup for reads
    with IUPAC ambiguity. Bowtie2 replaces those symbols with N; restoring them
    in the SAM strand keeps candidate sequences consistent with the read pool.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    if timeout is not None and (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
                                or not math.isfinite(timeout) or timeout <= 0):
        raise ValueError("timeout must be finite and positive, or None")
    if not read_paths:
        raise ValueError("At least one FASTQ path is required for alignment")
    if any("," in str(path) for path in read_paths):
        raise ValueError("Bowtie2 input paths cannot contain commas")
    strands = {None: None, "v": "--norc", "hv": "--norc", "lv": "--norc",
               "c": "--nofw", "hc": "--nofw", "lc": "--nofw",
               "--norc": "--norc", "--nofw": "--nofw"}
    if strand not in strands:
        raise ValueError("strand must be v/--norc or c/--nofw")
    preset = "--very-sensitive-local" if local else "--very-sensitive"
    command = [str(bowtie), "--no-unal", "--no-hd", preset, "--seed", "0",
               "--reorder", "--threads", str(threads)]
    if local:
        command.extend(["-L", "10", "-N", "1"])
    if strands[strand]:
        command.append(strands[strand])
    command.extend(["-x", str(index_prefix), "-U", ",".join(map(str, read_paths))])
    stderr_path = Path(stderr_path)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str], dict] = {}
    mapped_reads = 0
    with stderr_path.open("wb") as error_log:
        try:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error_log,
                                       text=True, encoding="ascii",
                                       start_new_session=(os.name == "posix"))
        except OSError as exc:
            raise AlignmentError(f"Cannot start Bowtie2: {exc}; log: {stderr_path}") from exc
        timed_out = threading.Event()
        finished = threading.Event()

        def expire():
            if not finished.is_set():
                timed_out.set()
                _stop_process_tree(process)

        timer = threading.Timer(timeout, expire) if timeout is not None else None
        if timer is not None:
            timer.daemon = True
            timer.start()
        try:
            assert process.stdout is not None
            for number, line in enumerate(process.stdout, 1):
                if line.startswith("@"):
                    if line[:3] not in {"@HD", "@SQ", "@RG", "@PG", "@CO"} or line[3:4] != "\t":
                        raise AlignmentError(f"Malformed SAM header at line {number}")
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 11:
                    raise AlignmentError(f"Malformed SAM at line {number}: expected at least 11 fields")
                try:
                    flag, position, mapq = int(fields[1]), int(fields[3]), int(fields[4])
                except ValueError as exc:
                    raise AlignmentError(f"Malformed SAM numeric field at line {number}") from exc
                if not (0 <= flag <= 65535 and 0 <= mapq <= 255 and position >= 0):
                    raise AlignmentError(f"Out-of-range SAM field at line {number}")
                if flag & (4 | 256 | 2048):
                    continue
                reference, sequence = fields[2], fields[9].upper().replace("U", "T")
                if (not fields[0] or reference in ("", "*") or position == 0 or not sequence
                        or not set(sequence) <= _IUPAC
                        or not re.fullmatch(r"(?:[1-9][0-9]*[MIDNSHP=X])+", fields[5])):
                    raise AlignmentError(f"Invalid mapped SAM record at line {number}")
                if fields[10] != "*" and (len(fields[10]) != len(sequence)
                        or _INVALID_QUALITY.search(fields[10])):
                    raise AlignmentError(f"Invalid SAM sequence/quality at line {number}")
                query_length = sum(int(length) for length, operation in
                                   re.findall(r"([0-9]+)([MIDNSHP=X])", fields[5])
                                   if operation in "MIS=X")
                if query_length != len(sequence):
                    raise AlignmentError(f"SAM CIGAR/sequence length mismatch at line {number}")
                if original_sequences:
                    original = original_sequences.get(fields[0])
                    if original is not None:
                        if not isinstance(original, str) or not set(original) <= _IUPAC:
                            raise AlignmentError(f"Invalid original IUPAC sequence for SAM read {fields[0]!r}")
                        restored = original.translate(_COMPLEMENT)[::-1] if flag & 16 else original
                        if len(restored) != len(sequence) or any(
                                emitted != expected and not (emitted == "N" and expected in _AMBIGUOUS)
                                for emitted, expected in zip(sequence, restored)):
                            raise AlignmentError(f"SAM/original sequence mismatch for read {fields[0]!r} at line {number}")
                        sequence = restored
                mapped_reads += 1
                key = sequence, reference
                if key not in grouped:
                    locus = re.search(r"(?:^|[|:_])(IGH|IGK|IGL|TRA|TRB|TRG|TRD)", reference)
                    grouped[key] = {"sequence": sequence, "count": 0, "reference": reference,
                                    "locus": locus.group(1) if locus else "", "mapq": mapq}
                grouped[key]["count"] += 1
                grouped[key]["mapq"] = max(grouped[key]["mapq"], mapq)
            returncode = process.wait()
        except BaseException as exc:
            cleanup_process(process)
            if timed_out.is_set() and isinstance(exc, Exception):
                detail = _stderr_tail(stderr_path)
                raise AlignmentError(f"Bowtie2 timed out after {timeout:g} seconds; log: {stderr_path}"
                                     + (f"\n{detail}" if detail else "")) from exc
            if isinstance(exc, UnicodeError):
                raise AlignmentError(f"Non-ASCII SAM from Bowtie2; log: {stderr_path}") from exc
            if isinstance(exc, AlignmentError):
                raise AlignmentError(f"{exc}; log: {stderr_path}") from exc
            raise
        finally:
            finished.set()
            if timer is not None:
                timer.cancel()
                timer.join()
            if process.stdout is not None:
                try:
                    process.stdout.close()
                except (OSError, ValueError):
                    pass
    if timed_out.is_set():
        detail = _stderr_tail(stderr_path)
        raise AlignmentError(f"Bowtie2 timed out after {timeout:g} seconds; log: {stderr_path}"
                             + (f"\n{detail}" if detail else ""))
    if returncode:
        detail = _stderr_tail(stderr_path)
        raise AlignmentError(f"Bowtie2 failed with exit {returncode}; log: {stderr_path}"
                             + (f"\n{detail}" if detail else ""))
    candidates = sorted(grouped.values(), key=lambda item:
                        (-item["count"], -item["mapq"], item["reference"], item["sequence"]))
    return MappingResult(candidates, mapped_reads, command)
