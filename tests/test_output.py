"""Real-aligner regressions against the two supplied human BCR examples."""
import json
import subprocess
import sys
import pytest

pytestmark = pytest.mark.integration


def read_chains(path):
    """Read all records, rejecting duplicate chain identifiers explicitly."""
    records, header, parts = {}, None, []

    def save():
        if header is not None:
            fields = header.split(";")
            assert len(fields) >= 2, f"Missing chain identifier: {header}"
            chain = fields[1]
            assert chain not in records, f"Duplicate chain record: {chain}"
            records[chain] = "".join(parts).upper()

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header, parts = line[1:], []
            else:
                assert header is not None, "Sequence precedes its header"
                parts.append(line)
    save()
    return records


def run_basic(repository_root, real_bowtie2, tmp_path, arguments):
    command = [sys.executable, "-m", "basic_assembler", "-b", real_bowtie2,
               "-d", str(repository_root / "db"), "-t", str(tmp_path),
               "-o", str(tmp_path / "out"), "-n", "result", "-p", "2", *arguments]
    return subprocess.run(command, cwd=repository_root, capture_output=True,
                          text=True, timeout=1200)


def test_se_human_bcr(repository_root, real_bowtie2, tmp_path):
    example = repository_root / "examples/cell_A1"
    result = run_basic(repository_root, real_bowtie2, tmp_path,
                       ["-SE", str(example / "ERR1421621_1M_reads.fastq.gz"), "-g", "human", "-i", "BCR"])
    assert result.returncode == 0, result.stdout + result.stderr
    actual = read_chains(tmp_path / "out/result.fasta")
    expected = read_chains(example / "result.fasta")
    assert set(actual) == set(expected) == {"heavy_chain", "light_chain"}
    assert actual == expected


def test_pe_human_bcr(repository_root, real_bowtie2, tmp_path):
    example = repository_root / "examples/H8_AW1"
    result = run_basic(repository_root, real_bowtie2, tmp_path,
                       ["-PE_1", str(example / "R1_100K_reads.fastq.gz"),
                        "-PE_2", str(example / "R2_100K_reads.fastq.gz"), "-g", "human", "-i", "BCR"])
    assert result.returncode == 0, result.stdout + result.stderr
    actual = read_chains(tmp_path / "out/result.fasta")
    sanger = read_chains(example / "sanger.fasta")
    assert set(actual) == set(sanger) == {"heavy_chain", "light_chain"}
    for chain, sequence in sanger.items():
        assert sequence in actual[chain]


@pytest.mark.parametrize("allow_partial", [False, True])
def test_empty_fastq_has_explicit_status(repository_root, real_bowtie2, tmp_path, allow_partial):
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    result = run_basic(repository_root, real_bowtie2, tmp_path,
                       ["-SE", str(empty), *(["-a"] if allow_partial else [])])
    if allow_partial:
        assert result.returncode == 0, result.stdout + result.stderr
        assert (tmp_path / "out/result.fasta").read_text() == ""
        report = json.loads((tmp_path / "out/result.run.json").read_text())
        assert report["status"] == "unsupported"
        assert report["read_statistics"]["consumed_reads"] == 0
    else:
        assert result.returncode != 0
        assert not (tmp_path / "out/result.fasta").exists()
        report = json.loads((tmp_path / "out/result.failed.json").read_text())
        assert report["status"] == "failed"
        assert report["error_type"] == "AssemblyIncomplete"
