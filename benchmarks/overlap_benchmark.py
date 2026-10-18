#!/usr/bin/env python3
"""Compare BASIC's scan-based overlap extension to its sorted-prefix index.

This measures one deterministic synthetic extension, not end-to-end assembly
or biological accuracy. It requires only Python and this checkout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import sys
import time

# Support ``python benchmarks/overlap_benchmark.py`` without installation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from basic_assembler.assembly import AssemblyOptions, OverlapIndex, _extend_paths


def legacy_scan_extend(seed, counts, read_length, options):
    """The original BASIC.py extension loop, with explicit safety bounds.

    Retains the full-pool substring scan, prefix scans, original entropy,
    shortest-overlap tie behavior, and insertion-order count ties. Changes
    are limited to empty-suffix handling and reporting/obeying the same step
    and length limits as the indexed comparison. The generated workload has
    only one supported path, so count ties cannot change the sequence.
    """
    sequence = seed
    steps = 0

    def result(termination, expanded):
        return {"sequence": sequence, "steps": steps,
                "expanded_states": expanded, "termination": termination}

    if len(sequence) > options.max_contig_length:
        raise ValueError("Seed exceeds max_contig_length")
    for expanded in range(1, options.max_steps + 1):
        suffixes = [sequence[-length:] for length in range(
            min(read_length, len(sequence)) - 1, options.min_overlap - 1, -1)]
        if not suffixes:
            return result("short_seed", expanded - 1)
        cache = {read: count for read, count in counts.items() if suffixes[-1] in read}
        entropies = [float("inf")]
        best_read = ""
        best_cut = 0
        for suffix in suffixes:
            matches = {read: count for read, count in cache.items() if read.startswith(suffix)}
            if not matches:
                continue
            denominator = float(sum(matches.values()) + 0.5)
            values = list(matches.values())
            probabilities = [float(value) / denominator for value in values]
            entropy = sum(-probability * math.log(probability, 2) for probability in probabilities)
            if entropy <= min(entropies) and max(values) >= options.min_support:
                best_read = list(matches)[values.index(max(values))]
                best_cut = -(len(best_read) - len(suffix))
                entropies.append(entropy)
        if len(entropies) == 1:
            return result("no_supported_overlap", expanded)
        if best_read.strip("N") in sequence:
            return result("repeated_read", expanded)
        candidate = (sequence + best_read[best_cut:]).strip("N")
        if len(candidate) > options.max_contig_length:
            return result("max_contig_length", expanded)
        sequence = candidate
        steps += 1
    return result("max_steps", options.max_steps)


def generate_workload(*, random_seed=7281, background_reads=100000,
                      read_length=100, contig_length=2500, min_overlap=15,
                      target_support=3):
    """Generate a reproducible pool with a unique supported target path."""
    if min_overlap < 1 or read_length <= min_overlap:
        raise ValueError("read_length must be greater than a positive min_overlap")
    if contig_length < read_length:
        raise ValueError("contig_length must be at least read_length")
    if background_reads < 0 or target_support < 2:
        raise ValueError("background_reads must be nonnegative and target_support at least 2")
    if contig_length - min_overlap + 1 > 4 ** min_overlap:
        raise ValueError("Target is too long to have unique minimum-length overlaps")
    target_read_count = contig_length - read_length + 1
    if background_reads + target_read_count > 4 ** read_length:
        raise ValueError("Requested more distinct reads than the sequence space permits")
    generator = random.Random(random_seed)
    for _ in range(100):
        reference = "".join(generator.choices("ACGT", k=contig_length))
        kmers = {reference[i:i + min_overlap] for i in range(contig_length - min_overlap + 1)}
        if len(kmers) == contig_length - min_overlap + 1:
            break
    else:
        raise ValueError("Could not generate a unique-overlap target; increase min_overlap")
    target_reads = {reference[i:i + read_length]: target_support
                    for i in range(target_read_count)}
    counts = {}
    attempts = 0
    while len(counts) < background_reads:
        attempts += 1
        if attempts > max(1000, background_reads * 50):
            raise ValueError("Could not generate enough distinct background reads")
        read = "".join(generator.choices("ACGT", k=read_length))
        if read not in target_reads:
            counts[read] = 1
    counts.update(target_reads)
    return reference[:read_length], reference, counts


def run_benchmark(*, random_seed=7281, background_reads=100000,
                  read_length=100, contig_length=2500, min_overlap=15,
                  target_support=3, max_steps=10000, max_contig_length=30000):
    parameters = {"random_seed": random_seed, "background_reads": background_reads,
                  "read_length": read_length, "contig_length": contig_length,
                  "min_overlap": min_overlap, "target_support": target_support,
                  "max_steps": max_steps, "max_contig_length": max_contig_length}
    options = AssemblyOptions(min_overlap=min_overlap, max_steps=max_steps,
                              max_contig_length=max_contig_length)
    if max_contig_length < contig_length:
        raise ValueError("max_contig_length must permit the full benchmark target")
    seed, reference, counts = generate_workload(**{
        key: value for key, value in parameters.items()
        if key not in {"max_steps", "max_contig_length"}})
    started = time.perf_counter()
    baseline = legacy_scan_extend(seed, counts, read_length, options)
    baseline_seconds = time.perf_counter() - started
    started = time.perf_counter()
    index = OverlapIndex(counts)
    index_seconds = time.perf_counter() - started
    started = time.perf_counter()
    indexed = _extend_paths(seed, index, options, {}, {})[0]
    extension_seconds = time.perf_counter() - started
    exact_output = baseline["sequence"] == indexed.sequence
    exact_reference = indexed.sequence == reference
    if not exact_output:
        raise AssertionError("Scan-based and indexed extension returned different sequences")
    if not exact_reference:
        raise AssertionError("Extension did not reconstruct the full synthetic reference; check bounds")
    digest = hashlib.sha256()
    for read in index.sequences:
        digest.update(f"{read}\t{index.counts[read]}\n".encode("ascii"))
    total_seconds = index_seconds + extension_seconds
    return {
        "benchmark": "synthetic_exact_overlap_extension",
        "parameters": parameters,
        "environment": {"python_version": platform.python_version(),
                        "python_implementation": platform.python_implementation(),
                        "platform": platform.platform(), "machine": platform.machine()},
        "workload": {"unique_reads": len(counts),
                     "supported_target_reads": contig_length - read_length + 1,
                     "seed_sequence": seed, "counts_sha256": digest.hexdigest(),
                     "assembled_sequence": indexed.sequence,
                     "sequence_sha256": hashlib.sha256(indexed.sequence.encode("ascii")).hexdigest()},
        "exact_output": exact_output,
        "exact_reference": exact_reference,
        "baseline": {"seconds": baseline_seconds, "steps": baseline["steps"],
                     "expanded_states": baseline["expanded_states"],
                     "termination": baseline["termination"]},
        "indexed": {"index_build_seconds": index_seconds,
                    "extension_seconds": extension_seconds,
                    "total_seconds": total_seconds, "steps": indexed.steps,
                    "termination": indexed.termination},
        "speedup_extension_only": baseline_seconds / extension_seconds if extension_seconds else None,
        "speedup_including_index": baseline_seconds / total_seconds if total_seconds else None,
        "measurement": "one run; scan first; generation and hashing excluded; no memory measurement",
        "scope": "unambiguous synthetic extension only; not end-to-end performance or biological validation",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in (("random-seed", 7281), ("background-reads", 100000),
                          ("read-length", 100), ("contig-length", 2500),
                          ("min-overlap", 15), ("target-support", 3),
                          ("max-steps", 10000), ("max-contig-length", 30000)):
        parser.add_argument("--" + name, type=int, default=default,
                            help=f"Default: {default}")
    parser.add_argument("--output", type=Path, help="Write JSON to this file instead of stdout")
    arguments = vars(parser.parse_args(argv))
    output = arguments.pop("output")
    try:
        report = run_benchmark(**arguments)
        serialized = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if output is None:
            print(serialized, end="")
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(serialized, encoding="utf-8")
    except (ValueError, AssertionError, OSError) as error:
        parser.exit(2, f"overlap_benchmark: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
