#!/usr/bin/env python3
"""Compare independently curated truth to assemblies with matching IDs.

Uses only the standard library. Run ``python benchmarks/benchmark.py --help``.
"""

import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records, identifier, parts = {}, None, []

    def save():
        if identifier is not None:
            sequence = "".join(parts).upper()
            if not sequence:
                raise ValueError("Empty FASTA sequence: " + identifier)
            records[identifier] = sequence

    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                words = line[1:].split()
                if not words:
                    raise ValueError("Empty FASTA identifier on line %d" % number)
                identifier = words[0]
                if identifier in records:
                    raise ValueError("Duplicate FASTA identifier: " + identifier)
                parts = []
            else:
                if identifier is None:
                    raise ValueError("FASTA sequence precedes a header")
                if any(char.isspace() for char in line):
                    raise ValueError("Whitespace inside FASTA sequence on line %d" % number)
                parts.append(line)
        save()
    return records


def _fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def compare_fastas(expected_path, actual_path):
    """Compare matching first-token FASTA IDs; orientation is not changed."""
    expected, actual = read_fasta(expected_path), read_fasta(actual_path)
    common = expected.keys() & actual.keys()
    exact = sorted(key for key in common if expected[key] == actual[key])
    contained = sorted(key for key in common if expected[key] in actual[key])
    partial = sorted(key for key in common if actual[key] in expected[key] and key not in exact)
    mismatched = sorted(common - set(contained) - set(partial))
    return {
        "expected_count": len(expected), "actual_count": len(actual),
        "common_id_count": len(common), "exact_match_count": len(exact),
        "exact_recall": _fraction(len(exact), len(expected)),
        "exact_precision": _fraction(len(exact), len(actual)),
        "expected_contained_count": len(contained),
        "containment_recall": _fraction(len(contained), len(expected)),
        "strict_partial_count": len(partial),
        "exact_ids": exact, "expected_contained_ids": contained,
        "strict_partial_ids": partial, "mismatched_ids": mismatched,
        "missing_ids": sorted(expected.keys() - actual.keys()),
        "extra_ids": sorted(actual.keys() - expected.keys()),
    }


def _read_annotations(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        if not reader.fieldnames or "sequence_id" not in reader.fieldnames:
            raise ValueError("Annotation TSV requires a sequence_id column")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Annotation TSV has duplicate column names")
        records = {}
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Annotation TSV contains a malformed row")
            identifier = row["sequence_id"]
            if not identifier or identifier in records:
                raise ValueError("Annotation sequence_id must be nonempty and unique")
            records[identifier] = row
    return records


def compare_annotations(expected_path, actual_path):
    """Evaluate non-null truth per field; missing actual values count as misses.

    Null truth is unscored. Junction and CDR3 are evaluated independently;
    neither is derived from the other or from a sequence motif.
    """
    expected, actual = _read_annotations(expected_path), _read_annotations(actual_path)
    output = {}
    for field in ("junction", "junction_aa", "cdr3", "cdr3_aa", "productive",
                  "v_call", "d_call", "j_call", "c_call", "locus"):
        truth = {key: row[field] for key, row in expected.items() if row.get(field)}
        if not truth:
            output[field] = {"truth_count": 0, "reported_count": 0,
                             "exact_match_count": 0, "exact_recall": None,
                             "accuracy_when_reported": None}
            continue
        reported = {key: actual.get(key, {}).get(field, "") for key in truth}
        count = sum(bool(value) for value in reported.values())
        exact = sum(reported[key] == value for key, value in truth.items())
        output[field] = {"truth_count": len(truth), "reported_count": count,
                         "exact_match_count": exact,
                         "exact_recall": _fraction(exact, len(truth)),
                         "accuracy_when_reported": _fraction(exact, count)}
    return output


def _execute(command, timeout):
    # A timed-out assembler can own alignment/worker processes. Keep the
    # command in its own POSIX session so those descendants are stopped too.
    process = subprocess.Popen(command, stdout=sys.stderr, stderr=sys.stderr,
                               start_new_session=(os.name == "posix"))
    try:
        return process.wait(timeout=timeout)
    except BaseException:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.wait()
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-fasta", required=True)
    parser.add_argument("--actual-fasta", required=True)
    parser.add_argument("--expected-airr")
    parser.add_argument("--actual-airr")
    parser.add_argument("--output", help="Write JSON report here instead of stdout")
    parser.add_argument("--label", help="Dataset/parameter label stored in the report")
    parser.add_argument("--timeout", type=float, default=3600,
                        help="Timeout in seconds for --run (default: 3600)")
    parser.add_argument("--run", nargs=argparse.REMAINDER,
                        help="Optional command and arguments to run before comparison; must be last")
    args = parser.parse_args(argv)
    if bool(args.expected_airr) != bool(args.actual_airr):
        parser.error("--expected-airr and --actual-airr must be supplied together")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    report = {"label": args.label, "matching": "first-token FASTA/sequence_id; same orientation"}
    try:
        if args.run is not None:
            if not args.run:
                parser.error("--run requires a command")
            started = time.monotonic()
            # Child output goes to stderr so stdout remains valid JSON.
            returncode = _execute(args.run, args.timeout)
            report["execution"] = {"command": args.run, "returncode": returncode,
                                   "wall_seconds": time.monotonic() - started}
            if returncode:
                raise RuntimeError("Benchmark command failed with exit status %d" % returncode)
        report["sequences"] = compare_fastas(args.expected_fasta, args.actual_fasta)
        if args.expected_airr:
            report["annotations"] = compare_annotations(args.expected_airr, args.actual_airr)
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        parser.exit(2, "benchmark: %s\n" % error)
    data = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(data, encoding="utf-8")
    else:
        print(data, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
