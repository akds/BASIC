from __future__ import annotations

import json

import pytest

from benchmarks.overlap_benchmark import generate_workload, legacy_scan_extend, main, run_benchmark
from basic_assembler.assembly import AssemblyOptions


def test_generated_pool_is_reproducible_and_has_one_supported_path():
    parameters = dict(random_seed=18, background_reads=30, read_length=20,
                      contig_length=140, min_overlap=7)
    first = generate_workload(**parameters)
    assert first == generate_workload(**parameters)
    seed, reference, counts = first
    assert seed == reference[:20]
    assert len(counts) == 30 + 140 - 20 + 1
    kmers = [reference[i:i + 7] for i in range(len(reference) - 6)]
    assert len(kmers) == len(set(kmers))
    assert {read for read, count in counts.items() if count >= 2} == {
        reference[i:i + 20] for i in range(len(reference) - 19)}


def test_small_benchmark_asserts_full_reference_and_records_metadata():
    report = run_benchmark(random_seed=18, background_reads=50, read_length=20,
                           contig_length=140, min_overlap=7)
    assert report["exact_output"]
    assert report["exact_reference"]
    assert len(report["workload"]["assembled_sequence"]) == 140
    assert report["environment"]["python_version"]
    assert report["baseline"]["termination"] == report["indexed"]["termination"]
    assert report["baseline"]["seconds"] >= 0
    assert report["indexed"]["total_seconds"] >= report["indexed"]["extension_seconds"]
    assert len(report["workload"]["counts_sha256"]) == 64


def test_scan_baseline_has_safe_step_and_length_bounds():
    counts = {"ACGT": 2, "CGTT": 2}
    step = legacy_scan_extend("AACG", counts, 4, AssemblyOptions(min_overlap=3, max_steps=1))
    length = legacy_scan_extend("AACG", counts, 4, AssemblyOptions(min_overlap=3, max_contig_length=5))
    assert step["sequence"] == length["sequence"] == "AACGT"
    assert step["termination"] == "max_steps"
    assert length["termination"] == "max_contig_length"


def test_benchmark_fails_instead_of_claiming_exactness_when_budget_is_too_small():
    with pytest.raises(AssertionError, match="full synthetic reference"):
        run_benchmark(background_reads=0, read_length=20, contig_length=140,
                      min_overlap=7, max_steps=1)


def test_cli_writes_json_and_creates_destination_directory(tmp_path):
    output = tmp_path / "results" / "overlap.json"
    assert main(["--background-reads", "20", "--read-length", "20",
                 "--contig-length", "100", "--min-overlap", "7",
                 "--output", str(output)]) == 0
    assert json.loads(output.read_text())["exact_output"]


def test_impossible_unique_overlap_workload_is_rejected():
    with pytest.raises(ValueError, match="unique minimum-length overlaps"):
        generate_workload(background_reads=0, read_length=5, contig_length=20, min_overlap=1)
