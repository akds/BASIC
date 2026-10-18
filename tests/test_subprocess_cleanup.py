"""Caller regressions for denied process-group signals during timeouts."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest

from basic_assembler import annotation, processes, reference_builder


BENCHMARK = Path(__file__).resolve().parents[1] / "benchmarks" / "benchmark.py"
SPEC = importlib.util.spec_from_file_location("basic_cleanup_benchmark", BENCHMARK)
benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark)


@pytest.mark.skipif(os.name != "posix", reason="Process-group signal fallback is POSIX-specific")
@pytest.mark.parametrize("caller", ["igblast", "builder", "benchmark"])
def test_group_permission_failure_preserves_timeout_and_reaps_child(monkeypatch, caller):
    children = []
    signals = []
    real_popen = subprocess.Popen

    def start(*args, **kwargs):
        process = real_popen(*args, **kwargs)
        children.append(process)
        return process

    def deny_group(pid, signal):
        signals.append((pid, signal))
        raise PermissionError("Process-group signal denied")

    monkeypatch.setattr(subprocess, "Popen", start)
    monkeypatch.setattr(processes.os, "killpg", deny_group)
    command = [sys.executable, "-c", "import time; time.sleep(60)"]
    timeout = 0.2
    try:
        with pytest.raises(subprocess.TimeoutExpired) as error:
            if caller == "igblast":
                annotation._run_igblast(command, os.environ.copy(), timeout)
            elif caller == "builder":
                reference_builder._run_builder(command, timeout)
            else:
                benchmark._execute(command, timeout)
        assert error.value.cmd == command
        assert error.value.timeout == timeout
        assert len(children) == 1
        process = children[0]
        assert signals and signals[0][0] == process.pid
        assert process.returncode is not None, "Direct child was not reaped"
        for stream in (process.stdin, process.stdout, process.stderr):
            assert stream is None or stream.closed
    finally:
        # A failed assertion must not leave the regression's sleeping child alive.
        for process in children:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=5)


def test_benchmark_source_cli_imports_cleanup_from_outside_checkout(tmp_path):
    completed = subprocess.run(
        [sys.executable, "-I", "-S", str(BENCHMARK), "--help"],
        cwd=tmp_path, capture_output=True, text=True, timeout=10,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--expected-fasta" in completed.stdout
