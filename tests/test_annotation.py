import hashlib
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from basic_assembler import annotation
from basic_assembler.annotation import (
    AnnotationError, REQUIRED_FIELDS, annotate_fasta, read_airr, write_airr,
    validate_annotation_options,
    annotation_reference_fingerprints,
)


def make_igblast(tmp_path, code):
    executable = tmp_path / "igblast with spaces;unused"
    executable.write_text("#!%s\nimport sys, pathlib, time\n"
                          "if '-version' in sys.argv:\n"
                          "    print('igblastn: test-version')\n"
                          "    sys.exit(0)\n" % sys.executable + code)
    executable.chmod(0o755)
    return executable


def query_file(tmp_path):
    query = tmp_path / "query ; $(unused).fa"
    query.write_text(">cell1 description\nACGT\n")
    return query


def output_code(rows):
    fields = list(REQUIRED_FIELDS) + ["locus", "cdr3", "cdr3_aa", "c_call"]
    text = "\t".join(fields) + "\n"
    text += "".join("\t".join(row.get(field, "") for field in fields) + "\n" for row in rows)
    return "pathlib.Path(sys.argv[sys.argv.index('-out')+1]).write_text(%r)\n" % text


def test_unannotated_airr_uses_nulls_and_preserves_metadata(tmp_path):
    output = tmp_path / "result.tsv"
    assert write_airr([{"sequence_id": "cell1", "sequence": "ACGT", "cell_id": "cell",
                        "status": "partial", "ambiguous": True, "read_support": 0,
                        "productive": None}], output) == 1
    row = read_airr(output)[0]
    assert set(REQUIRED_FIELDS).issubset(row)
    assert row["basic_status"] == "partial"
    assert row["basic_ambiguous"] == "T"
    assert row["basic_read_support"] == "0"
    assert row["productive"] == row["junction"] == row["junction_aa"] == ""
    assert "consensus_count" not in row


def test_airr_rejects_unsafe_cells_and_duplicate_ids_atomically(tmp_path):
    output = tmp_path / "existing.tsv"
    output.write_text("original")
    for rows in ([{"sequence_id": "bad\tidentifier"}],
                 [{"sequence_id": "one", "productive": "false"}],
                 [{"sequence_id": "one", "v_call": ["IGHV1", "IGHV2"]}],
                 [{"sequence_id": "one"}, {"sequence_id": "one"}]):
        with pytest.raises(ValueError):
            write_airr(rows, output)
        assert output.read_text() == "original"


def test_annotation_safe_command_and_metadata_merge(tmp_path):
    query = query_file(tmp_path)
    code = output_code([{"sequence_id": "cell1", "sequence": "ACGT", "productive": "T",
                         "junction": "TGTACGTGG", "junction_aa": "CTW", "cdr3": "ACG",
                         "cdr3_aa": "T", "locus": "IGH", "c_call": "IGHG1*01"}])
    code += "assert sys.argv[sys.argv.index('-outfmt')+1] == '19'\n"
    code += "assert sys.argv[sys.argv.index('-ig_seqtype')+1] == 'Ig'\n"
    code += "assert sys.argv[sys.argv.index('-germline_db_V')+1] == 'V ; $(unused)'\n"
    code += "assert sys.argv[sys.argv.index('-c_region_db')+1] == 'C database'\n"
    executable = make_igblast(tmp_path, code)
    output = tmp_path / "annotated.tsv"
    metadata = annotate_fasta(query, output, executable=executable, germline_v="V ; $(unused)",
                              germline_j="J", germline_d="D", germline_c="C database",
                              records=[{"sequence_id": "cell1", "cell_id": "cell",
                                        "status": "complete", "productive": False,
                                        "sequence": "WRONG", "ambiguous": False,
                                        "locus": "IGK", "basic_locus_source": "seed_alignment"}])
    row = read_airr(output)[0]
    assert row["productive"] == "T"
    assert row["sequence"] == "ACGT"
    assert row["junction"] != row["cdr3"]
    assert row["cell_id"] == "cell" and row["basic_status"] == "complete"
    assert row["c_call"] == "IGHG1*01"
    assert row["locus"] == "IGH"
    assert row["basic_seed_locus"] == "IGK"
    assert row["basic_locus_source"] == "igblast"
    assert metadata["record_count"] == 1
    assert metadata["version"] == "igblastn: test-version"
    assert metadata["warnings"]
    assert not list(tmp_path.glob(".basic-igblast-*"))


@pytest.mark.parametrize("code,match", [
    ("sys.stderr.write('database missing');sys.exit(4)\n", "database missing"),
    ("time.sleep(2)\n", "timeout"),
    ("pass\n", "without producing"),
    ("pathlib.Path(sys.argv[sys.argv.index('-out')+1]).write_text('bad\\n')\n", "required columns"),
    (output_code([{"sequence_id": "wrong", "sequence": "ACGT"}]), "do not match"),
    (output_code([]), "do not match"),
])
def test_annotation_failures_keep_previous_output(tmp_path, code, match):
    output = tmp_path / "result.tsv"
    output.write_text("prior result")
    with pytest.raises(AnnotationError, match=match):
        annotate_fasta(query_file(tmp_path), output, executable=make_igblast(tmp_path, code),
                       germline_v="V", germline_j="J", timeout=0.1 if match == "timeout" else 10)
    assert output.read_text() == "prior result"
    assert not list(tmp_path.glob(".basic-igblast-*"))


def test_tcr_igdata_and_auxiliary_options(tmp_path, monkeypatch):
    igdata = tmp_path / "igblast data"
    (igdata / "internal_data").mkdir(parents=True)
    auxiliary = igdata / "matched.aux"
    auxiliary.write_text("annotation fixture")
    monkeypatch.delenv("IGDATA", raising=False)
    code = "import os\nassert os.environ['IGDATA'] == %r\n" % str(igdata)
    code += "assert sys.argv[sys.argv.index('-ig_seqtype')+1] == 'TCR'\n"
    code += "assert sys.argv[sys.argv.index('-auxiliary_data')+1] == %r\n" % str(auxiliary)
    code += output_code([{"sequence_id": "cell1", "sequence": "ACGT"}])
    metadata = annotate_fasta(query_file(tmp_path), tmp_path / "result.tsv",
                              executable=make_igblast(tmp_path, code), receptor="TCR",
                              germline_v="V", germline_j="J", germline_d="D", auxiliary_data="matched.aux",
                              igdata=igdata)
    assert not metadata["warnings"]
    assert "IGDATA" not in os.environ


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group cleanup")
def test_timeout_reaps_descendant_processes(tmp_path):
    marker = tmp_path / "orphan-was-not-killed"
    child = "import pathlib,time;time.sleep(1.5);pathlib.Path(%r).touch()" % str(marker)
    code = "import subprocess\nsubprocess.Popen([sys.executable, '-c', %r])\ntime.sleep(5)\n" % child
    with pytest.raises(AnnotationError, match="timeout"):
        annotate_fasta(query_file(tmp_path), tmp_path / "result.tsv",
                       executable=make_igblast(tmp_path, code), germline_v="V", germline_j="J",
                       timeout=0.75)
    time.sleep(1.5)
    assert not marker.exists()


def test_empty_query_writes_header_without_igblast(tmp_path):
    query = tmp_path / "empty.fa"
    query.write_text("\n   \n")
    output = tmp_path / "empty.tsv"
    metadata = annotate_fasta(query, output, executable="not-installed", germline_v="V", germline_j="J")
    assert metadata["status"] == "skipped_empty_fasta"
    assert read_airr(output) == []


def test_empty_query_rejects_nonempty_assembly_metadata(tmp_path):
    query = tmp_path / "empty.fa"
    query.write_text("")
    with pytest.raises(ValueError, match="absent from input"):
        annotate_fasta(query, tmp_path / "empty.tsv", germline_v="V", germline_j="J",
                       records=[{"sequence_id": "unexpected", "sequence": "ACGT"}])


def test_annotation_preflight_checks_optional_resources(tmp_path):
    executable = make_igblast(tmp_path, "pass\n")
    assert validate_annotation_options(executable=executable, germline_v="V", germline_j="J") == str(executable)
    with pytest.raises(AnnotationError, match="executable"):
        validate_annotation_options(executable=tmp_path / "missing", germline_v="V", germline_j="J")
    with pytest.raises(ValueError, match="internal_data"):
        validate_annotation_options(executable=executable, germline_v="V", germline_j="J", igdata=tmp_path)
    with pytest.raises(ValueError, match="auxiliary file"):
        validate_annotation_options(executable=executable, germline_v="V", germline_j="J",
                                    auxiliary_data=tmp_path / "missing.aux")


def test_annotation_preflight_does_not_require_d_or_auxiliary(tmp_path):
    executable = make_igblast(tmp_path, "pass\n")
    assert validate_annotation_options(executable=executable, germline_v="V", germline_j="J",
                                       germline_d=None, auxiliary_data=None, receptor="TCR")


def test_annotation_preflight_rejects_absent_explicit_database_prefix(tmp_path):
    executable = make_igblast(tmp_path, "pass\n")
    prefix = tmp_path / "germline"
    with pytest.raises(ValueError, match="V database prefix"):
        validate_annotation_options(executable=executable, germline_v=prefix, germline_j="J")
    # Preflight confirms discoverability, not structural database validity.
    prefix.with_suffix(".nin").write_bytes(b"fixture index")
    assert validate_annotation_options(executable=executable, germline_v=prefix, germline_j="J")


def test_annotation_reference_fingerprints_include_only_selected_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prefix = tmp_path / "germline.v"
    for suffix in ("nhr", "nin", "nsq", "00.nin"):
        (tmp_path / (prefix.name + "." + suffix)).write_bytes(suffix.encode())
    (tmp_path / "germline.v.notes").write_text("not a BLAST database file")
    human = tmp_path / "internal_data" / "human"
    mouse = tmp_path / "internal_data" / "mouse"
    human.mkdir(parents=True)
    mouse.mkdir()
    (human / "human.ndm.imgt").write_text("human delineation")
    (mouse / "mouse.ndm.imgt").write_text("mouse delineation")
    auxiliary = tmp_path / "human.aux"
    auxiliary.write_text("J frames")
    result = annotation_reference_fingerprints(germline_v=prefix, germline_j="missing_J",
                                               auxiliary_data=auxiliary, igdata=tmp_path)
    files = {Path(item["path"]).name: item for item in result["files"]}
    assert len(files) == 6
    assert "germline.v.notes" not in files and "mouse.ndm.imgt" not in files
    assert files["human.aux"]["sha256"] == hashlib.sha256(b"J frames").hexdigest()
    assert files["germline.v.nin"]["roles"] == ["germline_v"]
    assert files["human.ndm.imgt"]["roles"] == ["internal_data"]
    assert any(item["role"] == "germline_j" for item in result["unresolved"])


def test_annotation_provenance_resolves_blastdb_aliases_and_detects_cycles(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    database = tmp_path / "database directory"
    database.mkdir()
    monkeypatch.setenv("BLASTDB", str(database))
    (database / "combined.nal").write_text('TITLE mixture\nDBLIST "part one" nested\nSEQIDLIST subset.ids\n')
    (database / "part one.nsq").write_bytes(b"sequences")
    (database / "nested.nal").write_text("DBLIST combined missing\n")
    (database / "subset.ids").write_text("allele1\n")
    result = annotation_reference_fingerprints(germline_v="combined", germline_j="combined")
    files = {Path(item["path"]).name: item for item in result["files"]}
    assert set(files) == {"combined.nal", "nested.nal", "part one.nsq", "subset.ids"}
    assert all(item["roles"] == ["germline_j", "germline_v"] for item in files.values())
    assert any("Cyclic" in item["reason"] for item in result["unresolved"])
    assert any(item["requested"] == "missing" for item in result["unresolved"])


def test_annotation_metadata_records_provenance_without_guessing_unresolved_paths(tmp_path):
    prefix = tmp_path / "V"
    (tmp_path / "V.nsq").write_bytes(b"sequence database")
    metadata = annotate_fasta(query_file(tmp_path), tmp_path / "annotated.tsv",
                              executable=make_igblast(tmp_path, output_code([{"sequence_id": "cell1", "sequence": "ACGT"}])),
                              germline_v=prefix, germline_j="unknown_J")
    assert metadata["references"]["files"][0]["sha256"] == hashlib.sha256(b"sequence database").hexdigest()
    assert any(item["requested"] == "unknown_J" for item in metadata["references"]["unresolved"])


def test_duplicate_fasta_ids_are_rejected(tmp_path):
    query = tmp_path / "bad.fa"
    query.write_text(">same\nACGT\n>same\nTGCA\n")
    with pytest.raises(ValueError, match="unique"):
        annotate_fasta(query, tmp_path / "result.tsv", germline_v="V", germline_j="J")


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf")])
def test_invalid_timeouts_are_rejected(tmp_path, timeout):
    with pytest.raises(ValueError, match="finite and positive"):
        annotate_fasta(query_file(tmp_path), tmp_path / "result.tsv",
                       germline_v="V", germline_j="J", timeout=timeout)


@pytest.mark.parametrize("text,match", [
    ("\t".join(REQUIRED_FIELDS) + "\ncell1\n", "field count"),
    ("\t".join(REQUIRED_FIELDS) + "\tsequence_id\n", "duplicate"),
])
def test_malformed_airr_is_rejected(tmp_path, text, match):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(AnnotationError, match=match):
        read_airr(path)


def test_airr_publish_failure_keeps_previous_file_and_cleans_temporary(tmp_path, monkeypatch):
    output = tmp_path / "previous.tsv"
    output.write_text("previous result")
    original_replace = annotation.os.replace

    def fail_publication(source, destination):
        if Path(destination) == output:
            raise OSError("simulated disk failure")
        return original_replace(source, destination)

    monkeypatch.setattr(annotation.os, "replace", fail_publication)
    with pytest.raises(OSError, match="disk failure"):
        write_airr([{"sequence_id": "cell1", "sequence": "ACGT"}], output)
    assert output.read_text() == "previous result"
    assert not list(tmp_path.glob(".previous.tsv.*"))


@pytest.mark.parametrize("rows, message", [
    ([{"sequence_id": "cell1", "sequence": "ACGT"}] * 2, "unique"),
    ([{"sequence_id": "cell1", "sequence": "ACGT", "productive": "True"}], "AIRR boolean"),
    ([{"sequence_id": "", "sequence": "ACGT"}], "nonempty"),
])
def test_invalid_annotation_rows_do_not_replace_prior_output(tmp_path, rows, message):
    output = tmp_path / "prior.tsv"
    output.write_text("prior result")
    with pytest.raises(AnnotationError, match=message):
        annotate_fasta(query_file(tmp_path), output,
                       executable=make_igblast(tmp_path, output_code(rows)),
                       germline_v="V", germline_j="J")
    assert output.read_text() == "prior result"
    assert not list(tmp_path.glob(".basic-igblast-*"))


def test_malformed_reference_alias_is_reported_as_unresolved(tmp_path):
    prefix = tmp_path / "V"
    (tmp_path / "V.nal").write_text('DBLIST "unterminated\n')
    provenance = annotation_reference_fingerprints(germline_v=prefix, germline_j="unknown")
    assert any("Cannot resolve alias" in row["reason"] for row in provenance["unresolved"])
    assert any(Path(row["path"]).name == "V.nal" for row in provenance["files"])


@pytest.mark.parametrize("content, message", [
    (">empty\n>next\nACGT\n", "Empty FASTA sequence"),
    (">empty\n", "Empty FASTA sequence"),
    (">   \nACGT\n", "identifier must not be empty"),
    ("ACGT\n", "precedes its header"),
])
def test_malformed_query_fasta_is_rejected_before_annotation(tmp_path, content, message):
    query = tmp_path / "invalid.fa"
    query.write_text(content)
    output = tmp_path / "previous.tsv"
    output.write_text("previous result")
    with pytest.raises(ValueError, match=message):
        annotate_fasta(query, output, executable="not-installed", germline_v="V", germline_j="J")
    assert output.read_text() == "previous result"


@pytest.mark.parametrize("records, message", [
    ([{"sequence_id": "cell1"}, {"sequence_id": "cell1"}], "unique"),
    ([{"sequence": "ACGT"}], "nonempty"),
    ([{"sequence_id": "other-cell"}], "absent from input"),
])
def test_invalid_metadata_join_does_not_start_annotation(tmp_path, records, message):
    marker = tmp_path / "annotation-started"
    code = "pathlib.Path(%r).touch()\n" % str(marker)
    executable = make_igblast(tmp_path, code)
    with pytest.raises(ValueError, match=message):
        annotate_fasta(query_file(tmp_path), tmp_path / "result.tsv", executable=executable,
                       germline_v="V", germline_j="J", records=records)
    assert not marker.exists()


def test_annotation_cannot_overwrite_query_via_symlink(tmp_path):
    query = query_file(tmp_path)
    output_alias = tmp_path / "query-link.tsv"
    output_alias.symlink_to(query)
    before = query.read_bytes()
    with pytest.raises(ValueError, match="must differ"):
        annotate_fasta(query, output_alias, germline_v="V", germline_j="J")
    assert query.read_bytes() == before
    assert output_alias.is_symlink()


@pytest.mark.parametrize("failure", ["unsupported", "unavailable", "timeout"])
def test_version_probe_failure_does_not_discard_valid_annotation(tmp_path, monkeypatch, failure):
    executable = make_igblast(tmp_path, output_code([{"sequence_id": "cell1", "sequence": "ACGT"}]))
    original_run = annotation._run_igblast

    def version_failure(command, environment, timeout):
        if "-version" not in command:
            return original_run(command, environment, timeout)
        if failure == "unavailable":
            raise OSError("executable disappeared after annotation")
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, timeout)
        return subprocess.CompletedProcess(command, 2, "unsupported option", "")

    monkeypatch.setattr(annotation, "_run_igblast", version_failure)
    output = tmp_path / "result.tsv"
    metadata = annotate_fasta(query_file(tmp_path), output, executable=executable,
                              germline_v="V", germline_j="J")
    assert metadata["status"] == "complete"
    assert metadata["version"] is None
    assert read_airr(output)[0]["sequence"] == "ACGT"


def test_relative_alias_dependencies_are_deduplicated_and_missing_filters_reported(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "aliases").mkdir()
    (tmp_path / "shared").mkdir()
    (tmp_path / "shared" / "V.nsq").write_bytes(b"sequences")
    (tmp_path / "aliases" / "combined.nal").write_text(
        '# relative dependencies\n\nDBLIST ../shared/V ../shared/V\nGILIST absent.ids\n')
    result = annotation_reference_fingerprints(germline_v="aliases/combined", germline_j="unknown")
    files = [item["path"] for item in result["files"]]
    assert files.count(str(tmp_path / "shared" / "V.nsq")) == 1
    assert len([item for item in result["database_prefixes"] if item["resolved_prefix"].endswith("/shared/V")]) == 1
    assert any(item["requested"].endswith("absent.ids") for item in result["unresolved"])


@pytest.mark.parametrize("option, value, message", [
    ("organism", "", "organism"),
    ("receptor", "unknown", "BCR or TCR"),
    ("germline_v", None, "V and J"),
    ("threads", 0, "positive integer"),
])
def test_preflight_rejects_invalid_api_options_before_executable_lookup(option, value, message):
    options = {"germline_v": "V", "germline_j": "J", option: value}
    with pytest.raises(ValueError, match=message):
        validate_annotation_options(executable="not-installed", **options)


def test_unlaunchable_executable_preserves_previous_annotation(tmp_path):
    executable = tmp_path / "broken-igblast"
    executable.write_text("#!/nonexistent/basic-test-interpreter\n")
    executable.chmod(0o755)
    output = tmp_path / "result.tsv"
    output.write_text("previous result")
    with pytest.raises(AnnotationError, match="Could not run IgBLAST"):
        annotate_fasta(query_file(tmp_path), output, executable=executable,
                       germline_v="V", germline_j="J")
    assert output.read_text() == "previous result"
    assert not list(tmp_path.glob(".basic-igblast-*"))
