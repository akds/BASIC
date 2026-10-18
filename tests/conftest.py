"""Shared fixtures keep external-tool integration explicit and optional."""
import os
from pathlib import Path
import shutil
import pytest
from basic_assembler.alignment import AlignmentError, resolve_bowtie


@pytest.fixture(scope="session")
def repository_root():
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def real_bowtie2():
    requested = os.environ.get("BASIC_BOWTIE2")
    executable = requested or shutil.which("bowtie2")
    if not executable:
        pytest.skip("Real Bowtie2 is unavailable; set BASIC_BOWTIE2 to enable integration tests")
    try:
        return resolve_bowtie(executable)[0]
    except AlignmentError as exc:
        if requested:
            pytest.fail(f"Configured BASIC_BOWTIE2 is unusable: {exc}")
        pytest.skip(f"Bowtie2 on PATH is unusable: {exc}")
