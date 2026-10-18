"""CLI configuration checks independent of a real aligner installation."""
import pytest
from basic_assembler.cli import config_args, parse_args


@pytest.fixture
def configured_input(tmp_path, repository_root, monkeypatch):
    fastq = tmp_path / "reads with spaces.fastq"
    fastq.write_text("@read\nACGTACGTTGCACTGA\n+\nIIIIIIIIIIIIIIII\n")
    monkeypatch.setattr("basic_assembler.alignment.resolve_bowtie",
                        lambda requested: ("/stub/bowtie2", "stub version"))
    return {"FASTQ": str(fastq), "database_path": str(repository_root / "db"),
            "output_location": str(tmp_path / "output not created")}


def test_single_and_paired_configuration(configured_input):
    single = config_args(configured_input)
    assert single["read_type"] == "single"
    assert single["single"] == [configured_input["FASTQ"]]
    paired = config_args(dict(configured_input, FASTQ="", LEFT=configured_input["FASTQ"],
                              RIGHT=configured_input["FASTQ"]))
    assert paired["read_type"] == "paired"
    assert paired["left"] == paired["right"]


@pytest.mark.parametrize("changes", [
    {"FASTQ": ""}, {"FASTQ": "", "LEFT": "missing"},
    {"FASTQ": "", "RIGHT": "missing"}, {"LEFT": "conflicting"},
    {"type": "invalid"}, {"name": "../escape"}, {"genome": "../escape"},
    {"num_threads": 0}, {"subsample": -1}, {"subsample": 0},
    {"trim_quality": 94}, {"min_support": 0},
    {"quality_weighted": True}, {"mate_aware": True, "scoring": "evidence"},
])
def test_invalid_options_fail_before_running(configured_input, changes):
    with pytest.raises(ValueError):
        config_args(dict(configured_input, **changes))


def test_validation_does_not_create_output_directory(configured_input):
    from pathlib import Path
    config_args(configured_input)
    assert not Path(configured_input["output_location"]).exists()


def test_missing_database_shard_is_rejected(configured_input, tmp_path):
    from basic_assembler.alignment import AlignmentError
    database = tmp_path / "incomplete" / "human" / "BCR"
    database.mkdir(parents=True)
    (database / "hc").write_text(">reference\nACGT\n")
    with pytest.raises(AlignmentError, match="Incomplete"):
        config_args(dict(configured_input, database_path=str(database.parent.parent)))


def test_custom_species_uses_database_directory(configured_input, tmp_path):
    database = tmp_path / "custom" / "test_species" / "BCR"
    database.mkdir(parents=True)
    for region in ("hv", "hc", "lv", "lc"):
        for shard in ("1", "2", "3", "4", "rev.1", "rev.2"):
            (database / f"{region}.{shard}.bt2").write_bytes(b"stub index")
    result = config_args(dict(configured_input, genome="test_species",
                              database_path=str(database.parent.parent)))
    assert result["db_path"] == str(database)


def test_unequal_lane_lists_are_rejected(configured_input):
    path = configured_input["FASTQ"]
    with pytest.raises(ValueError, match="same number"):
        config_args(dict(configured_input, FASTQ="", LEFT=path, RIGHT=f"{path},{path}"))


def test_default_algorithm_options_are_explicit():
    args = parse_args([])
    assert args["scoring"] == "legacy"
    assert args["beam_width"] == args["max_anchors"] == args["max_candidates"] == 1
    assert args["subsample"] is None
    assert not args["quality_weighted"]
    assert not args["mate_aware"]
