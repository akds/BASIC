import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "benchmark.py"
SPEC = importlib.util.spec_from_file_location("basic_benchmark", SCRIPT)
benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark)


def test_sequence_metrics_distinguish_full_partial_and_mismatch(tmp_path):
    truth = tmp_path / "truth.fa"
    actual = tmp_path / "actual.fa"
    truth.write_text(">exact desc\nACGT\n>extended\nATCG\n>partial\nGGCA\n>wrong\nGATT\n>missing\nTCTA\n")
    actual.write_text(">exact\nacgt\n>extended\nAATCGA\n>partial\nGC\n>wrong\nAATT\n>extra\nGGGG\n")
    metrics = benchmark.compare_fastas(truth, actual)
    assert metrics["exact_match_count"] == 1
    assert metrics["exact_recall"] == metrics["exact_precision"] == 0.2
    assert metrics["expected_contained_count"] == 2
    assert metrics["strict_partial_ids"] == ["partial"]
    assert metrics["mismatched_ids"] == ["wrong"]
    assert metrics["missing_ids"] == ["missing"]
    assert metrics["extra_ids"] == ["extra"]


def test_annotation_metrics_score_null_truth_and_predictions_correctly(tmp_path):
    truth = tmp_path / "truth.tsv"
    actual = tmp_path / "actual.tsv"
    truth.write_text("sequence_id\tjunction\tcdr3\tproductive\n"
                     "one\tTGTACGTGG\tACG\tT\ntwo\tTGTAAATGG\tAAA\tF\nunknown\t\t\t\n")
    actual.write_text("sequence_id\tjunction\tcdr3\tproductive\n"
                      "one\tTGTACGTGG\tTGTACGTGG\tT\ntwo\t\t\t\nunknown\tNNN\tNNN\tF\n")
    metrics = benchmark.compare_annotations(truth, actual)
    assert metrics["junction"]["exact_recall"] == 0.5
    assert metrics["junction"]["accuracy_when_reported"] == 1
    assert metrics["cdr3"]["exact_recall"] == 0
    assert metrics["productive"]["truth_count"] == 2
    assert metrics["v_call"]["exact_recall"] is None


def test_empty_truth_is_not_perfect_accuracy(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    metrics = benchmark.compare_fastas(empty, empty)
    assert metrics["exact_recall"] is None
    assert metrics["exact_precision"] is None


@pytest.mark.parametrize("content", [">duplicate\nAC\n>duplicate\nGT\n", ">empty\n", "ACGT\n"])
def test_invalid_fasta_is_rejected(tmp_path, content):
    fasta = tmp_path / "bad.fa"
    fasta.write_text(content)
    with pytest.raises(ValueError):
        benchmark.read_fasta(fasta)


def test_cli_executes_and_times_command_without_polluting_json(tmp_path):
    truth = tmp_path / "truth.fa"
    actual = tmp_path / "actual.fa"
    truth.write_text(">id\nACGT\n")
    command = [sys.executable, str(SCRIPT), "--expected-fasta", str(truth),
               "--actual-fasta", str(actual), "--run", sys.executable, "-c",
               "import pathlib; pathlib.Path(%r).write_text('>id\\nACGT\\n');print('child output')" % str(actual)]
    completed = subprocess.run(command, capture_output=True, text=True, check=True)
    result = json.loads(completed.stdout)
    assert "child output" in completed.stderr
    assert result["sequences"]["exact_recall"] == 1
    assert result["execution"]["wall_seconds"] > 0


def test_cli_failure_does_not_score_stale_outputs(tmp_path):
    fasta = tmp_path / "existing.fa"
    fasta.write_text(">id\nACGT\n")
    completed = subprocess.run([sys.executable, str(SCRIPT), "--expected-fasta", str(fasta),
                                "--actual-fasta", str(fasta), "--run", sys.executable,
                                "-c", "raise SystemExit(7)"], capture_output=True, text=True)
    assert completed.returncode == 2
    assert "exit status 7" in completed.stderr
    assert completed.stdout == ""
