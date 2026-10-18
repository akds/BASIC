"""Synthetic fixtures exercise software behavior, not empirical accuracy."""
import json
import subprocess
import sys
import pytest
from basic_assembler.reads import load_read_pool
from benchmarks.benchmark import read_fasta
from benchmarks.generate_synthetic import _references, generate_cell


@pytest.mark.parametrize("species,receptor", [("human", "BCR"), ("human", "TCR"),
                                               ("mouse", "BCR"), ("mouse", "TCR")])
def test_reference_constructs_have_explicit_truth(tmp_path, species, receptor):
    manifest = generate_cell(tmp_path, species=species, receptor=receptor)
    assert manifest["dataset_kind"] == "synthetic_software_validation"
    truth = read_fasta(tmp_path / "truth.fasta")
    assert len(truth) == 2
    assert truth == {item["sequence_id"]: item["sequence"] for item in manifest["constructs"]}
    for item in manifest["constructs"]:
        variable = manifest["references"][item["variable_region"]]
        constant = manifest["references"][item["constant_region"]]
        assert variable["locus"] == constant["locus"] == item["locus"]
        assert constant["construct_orientation"] == "reverse_complement"
        assert variable["record_index"] == constant["record_index"] == 0
        assert len(variable["source_sha256"]) == 64
        assert item["sequence"][item["junction_start"]:item["junction_end"]] == item["synthetic_junction"]
    if species == "mouse" and receptor == "TCR":
        assert [item["chain"] for item in manifest["constructs"]] == ["TRB", "TRA"]


def test_generator_is_reproducible_with_variable_lengths_noise_and_multiple_candidates(tmp_path):
    options = {"read_lengths": (50, 75), "copies": 1, "step": 40,
               "substitution_rate": 0.03, "noise_reads": 7, "seed": 42,
               "candidates_per_chain": 2}
    first = generate_cell(tmp_path / "one", **options)
    second = generate_cell(tmp_path / "two", **options)
    assert first == second
    for name in ("reads.fastq", "truth.fasta", "manifest.json", "samples.tsv"):
        assert (tmp_path / "one" / name).read_bytes() == (tmp_path / "two" / name).read_bytes()
    assert len(first["constructs"]) == 4
    lines = (tmp_path / "one/reads.fastq").read_text().splitlines()
    assert {len(sequence) for sequence in lines[1::4]} == {50, 75}
    assert sum(header.startswith("@synthetic_noise") for header in lines[::4]) == 7
    assert any("+" in quality for quality in lines[3::4])


def test_paired_fixture_is_synchronized(tmp_path):
    manifest = generate_cell(tmp_path, paired=True, read_lengths=(50, 75), copies=2)
    pool = load_read_pool([], [str(tmp_path / "R1.fastq")], [str(tmp_path / "R2.fastq")], retain_mates=True)
    assert pool.stats["consumed_reads"] == manifest["physical_reads"]
    assert pool.stats["consumed_pairs"] * 2 == manifest["physical_reads"]
    assert pool.mate_links


def test_reference_choice_is_checked_and_cannot_overwrite_fixture(tmp_path):
    with pytest.raises(ValueError, match="reference index"):
        generate_cell(tmp_path / "bad", reference_indices={"hv": -1})
    generate_cell(tmp_path / "good")
    with pytest.raises(FileExistsError):
        generate_cell(tmp_path / "good")


def gamma_delta_indices(repository_root):
    choices = {}
    for region, locus in (("hv", "TRD"), ("hc", "TRD"), ("lv", "TRG"), ("lc", "TRG")):
        for index, (header, sequence) in enumerate(_references(repository_root / "db/mouse/TCR" / region)):
            if f"|{locus}" in header and sequence.strip("N") and set(sequence.strip("N")) <= set("ACGT"):
                choices[region] = index
                break
        assert region in choices, f"No unambiguous {locus} reference in {region}"
    return choices


def test_mouse_gamma_delta_selection_is_explicit(tmp_path, repository_root):
    choices = gamma_delta_indices(repository_root)
    manifest = generate_cell(tmp_path, species="mouse", receptor="TCR", reference_indices=choices)
    assert [item["chain"] for item in manifest["constructs"]] == ["TRD", "TRG"]
    assert {region: value["record_index"] for region, value in manifest["references"].items()} == choices


def assemble_fixture(tmp_path, repository_root, real_bowtie2, manifest, allow_partial=False):
    dataset = tmp_path / "synthetic"
    output = tmp_path / "output"
    reads = (["-PE_1", str(dataset / "R1.fastq"), "-PE_2", str(dataset / "R2.fastq")]
             if manifest["options"]["paired"] else ["-SE", str(dataset / "reads.fastq")])
    command = [sys.executable, "-m", "basic_assembler", *reads,
               "-g", manifest["species"], "-i", manifest["receptor"], "-d", str(repository_root / "db"),
               "-b", real_bowtie2, "-n", manifest["sample_id"], "-o", str(output), "-p", "2"]
    if allow_partial:
        command.append("-a")
    result = subprocess.run(command, cwd=repository_root, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    actual = read_fasta(output / f"{manifest['sample_id']}.fasta")
    report = json.loads((output / f"{manifest['sample_id']}.run.json").read_text())
    return actual, report


@pytest.mark.integration
@pytest.mark.parametrize("species,receptor,gamma_delta", [
    ("human", "BCR", False), ("human", "TCR", False),
    ("mouse", "BCR", False), ("mouse", "TCR", False),
    pytest.param("mouse", "TCR", True, id="mouse-gamma-delta"),
])
def test_clean_synthetic_assembly_matches_constructs(tmp_path, repository_root, real_bowtie2, species, receptor, gamma_delta):
    choices = gamma_delta_indices(repository_root) if gamma_delta else None
    manifest = generate_cell(tmp_path / "synthetic", species=species, receptor=receptor,
                             step=15, copies=3, reference_indices=choices)
    actual, report = assemble_fixture(tmp_path, repository_root, real_bowtie2, manifest)
    assert actual == read_fasta(tmp_path / "synthetic/truth.fasta")
    assert report["status"] == "complete"


@pytest.mark.integration
def test_clean_variable_length_paired_reads_match_constructs(tmp_path, repository_root, real_bowtie2):
    manifest = generate_cell(tmp_path / "synthetic", read_lengths=(50, 75), paired=True, step=15, copies=3)
    actual, report = assemble_fixture(tmp_path, repository_root, real_bowtie2, manifest)
    assert actual == read_fasta(tmp_path / "synthetic/truth.fasta")
    assert report["read_statistics"]["consumed_pairs"] * 2 == manifest["physical_reads"]


@pytest.mark.integration
def test_insufficient_synthetic_support_is_reported_as_partial(tmp_path, repository_root, real_bowtie2):
    manifest = generate_cell(tmp_path / "synthetic", copies=1, step=30)
    actual, report = assemble_fixture(tmp_path, repository_root, real_bowtie2, manifest, allow_partial=True)
    truth = read_fasta(tmp_path / "synthetic/truth.fasta")
    assert report["status"] == "partial"
    assert actual
    assert all(any(sequence in full for full in truth.values()) for sequence in actual.values())
