from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.parse import parse_qs, urlsplit
from urllib.request import Request
from urllib.error import URLError

import pytest

from basic_assembler import reference_builder as builder
from basic_assembler.assembly import revcom


SEQUENCE = "ACGTTAGCTAGGCTAACGAT"


@pytest.fixture
def fake_builder(tmp_path):
    executable = tmp_path / "builder with spaces"
    executable.write_text(f"#!{sys.executable}\n" + """
import os
from pathlib import Path
import sys
if '--version' in sys.argv:
    print('test bowtie2-build 1.0')
    raise SystemExit(0)
prefix = Path(sys.argv[2])
if os.environ.get('BASIC_TEST_SLEEP') == prefix.name:
    import time
    time.sleep(2)
if os.environ.get('BASIC_TEST_BUILD_FAIL') == prefix.name:
    print('deliberate build failure')
    raise SystemExit(9)
suffixes = ['1', '2', '3', '4', 'rev.1', 'rev.2']
if os.environ.get('BASIC_TEST_OMIT_SHARD') == prefix.name:
    suffixes.pop()
for suffix in suffixes:
    Path(str(prefix) + '.' + suffix + '.bt2').write_bytes(b'index')
""")
    executable.chmod(0o755)
    return executable


@pytest.fixture
def manifest(tmp_path):
    regions = {}
    for region in builder.REGIONS:
        (tmp_path / f"{region}.fa").write_text(f">gene_{region} reference description\n{SEQUENCE}\n")
        regions[region] = [{"path": f"{region}.fa", "orientation": "forward",
                            "locus": "IGH" if region.startswith("h") else "IGK"}]
    path = tmp_path / "sources.json"
    path.write_text(json.dumps({"schema_version": 1, "regions": regions,
                                "metadata": {"scientific_name": "Example species", "taxon_id": 123}}))
    return path


def update_manifest(path, update):
    data = json.loads(path.read_text())
    update(data)
    path.write_text(json.dumps(data))


def build(manifest, fake_builder, tmp_path, **kwargs):
    return builder.build_references(manifest, tmp_path / "references", "example_species", "BCR",
                                     fake_builder, **kwargs)


def snapshot(path):
    return {str(file.relative_to(path)): file.read_bytes() for file in path.rglob("*") if file.is_file()}


class FakeResponse(io.BytesIO):
    def __init__(self, data, url, headers=None):
        super().__init__(data)
        self.url = url
        self.headers = headers or {}

    def geturl(self):
        return self.url


class FakeOpener:
    def __init__(self, contents):
        self.contents = contents
        self.requests = []

    def open(self, request, timeout):
        self.requests.append((request.full_url, timeout))
        return FakeResponse(self.contents, request.full_url)


def test_build_validates_indexes_and_orients_constants(manifest, fake_builder, tmp_path):
    result = build(manifest, fake_builder, tmp_path)
    target = Path(result["path"])
    assert (target / "hv").read_text().splitlines() == [">IGH|gene_hv reference description", SEQUENCE]
    assert (target / "hc").read_text().splitlines()[1] == revcom(SEQUENCE)
    assert (target / "lc").read_text().splitlines()[1] == revcom(SEQUENCE)
    assert len(list(target.glob("*.bt2"))) == 24
    recorded = json.loads((target / "reference_manifest.json").read_text())
    assert recorded["builder"]["version"] == "test bowtie2-build 1.0"
    assert recorded["metadata"]["taxon_id"] == 123
    source = recorded["regions"]["hc"]["sources"][0]
    assert source["sha256"] == hashlib.sha256((manifest.parent / "hc.fa").read_bytes()).hexdigest()
    assert source["reverse_complemented"]
    assert source["orientation"] == "forward"
    assert not list(target.parent.glob(".*building-*"))
    assert not list(target.parent.glob("*.lock"))


def test_as_is_constant_is_not_reverse_complemented(manifest, fake_builder, tmp_path):
    update_manifest(manifest, lambda data: data["regions"]["hc"][0].update(orientation="as_is"))
    result = build(manifest, fake_builder, tmp_path)
    assert (Path(result["path"]) / "hc").read_text().splitlines()[1] == SEQUENCE


@pytest.mark.parametrize("failure_variable", ["BASIC_TEST_BUILD_FAIL", "BASIC_TEST_OMIT_SHARD"])
def test_failed_overwrite_preserves_existing_reference(manifest, fake_builder, tmp_path, monkeypatch, failure_variable):
    target = Path(build(manifest, fake_builder, tmp_path)["path"])
    before = snapshot(target)
    (manifest.parent / "hv.fa").write_text(">replacement\nATGCGTACGTAGCTAGCGTA\n")
    monkeypatch.setenv(failure_variable, "lc")
    with pytest.raises(builder.ReferenceBuildError):
        build(manifest, fake_builder, tmp_path, overwrite=True)
    assert snapshot(target) == before
    assert not list(target.parent.glob(".*building-*"))


def test_failed_install_restores_old_directory(manifest, fake_builder, tmp_path, monkeypatch):
    target = Path(build(manifest, fake_builder, tmp_path)["path"])
    before = snapshot(target)
    original_replace = builder.os.replace

    def fail_new_directory(source, destination):
        if ".building-" in Path(source).name:
            raise OSError("simulated installation failure")
        return original_replace(source, destination)

    monkeypatch.setattr(builder.os, "replace", fail_new_directory)
    with pytest.raises(builder.ReferenceBuildError, match="installation failure"):
        build(manifest, fake_builder, tmp_path, overwrite=True)
    assert snapshot(target) == before


def test_build_timeout_preserves_existing_reference(manifest, fake_builder, tmp_path, monkeypatch):
    target = Path(build(manifest, fake_builder, tmp_path)["path"])
    before = snapshot(target)
    monkeypatch.setenv("BASIC_TEST_SLEEP", "hc")
    with pytest.raises(builder.ReferenceBuildError, match="timed out"):
        build(manifest, fake_builder, tmp_path, overwrite=True, build_timeout=0.1)
    assert snapshot(target) == before


def test_overwrite_is_explicit_and_successful_replacement_is_complete(manifest, fake_builder, tmp_path):
    target = Path(build(manifest, fake_builder, tmp_path)["path"])
    with pytest.raises(builder.ReferenceBuildError, match="--overwrite"):
        build(manifest, fake_builder, tmp_path)
    (manifest.parent / "hv.fa").write_text(">replacement\nATGCGTACGTAGCTAGCGTA\n")
    build(manifest, fake_builder, tmp_path, overwrite=True)
    assert "replacement" in (target / "hv").read_text()
    assert len(list(target.glob("*.bt2"))) == 24
    assert not list(target.parent.glob(".*backup-*"))


def test_https_download_hash_and_provenance(manifest, fake_builder, tmp_path, monkeypatch):
    content = b">remote_gene\nACGTTAGCTAGGCTAACGAT\n"
    opener = FakeOpener(content)
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)

    def use_remote(data):
        source = data["regions"]["hv"][0]
        source.pop("path")
        source.update(url="https://references.example/variable.fa",
                      sha256=hashlib.sha256(content).hexdigest())

    update_manifest(manifest, use_remote)
    result = build(manifest, fake_builder, tmp_path)
    source = result["manifest"]["regions"]["hv"]["sources"][0]
    assert source["resolved_url"] == "https://references.example/variable.fa"
    assert source["sha256"] == hashlib.sha256(content).hexdigest()
    assert opener.requests == [(source["url"], 30.0)]


def test_checksum_mismatch_fails_before_install(manifest, fake_builder, tmp_path):
    update_manifest(manifest, lambda data: data["regions"]["hv"][0].update(sha256="0" * 64))
    with pytest.raises(builder.ReferenceBuildError, match="SHA256 mismatch"):
        build(manifest, fake_builder, tmp_path)
    assert not (tmp_path / "references/example_species/BCR").exists()


def test_accession_uses_ncbi_nucleotide_fasta_endpoint(manifest, fake_builder, tmp_path, monkeypatch):
    opener = FakeOpener(b">NM_123456.2 example\nACGTTAGCTAGGCTAACGAT\n")
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)

    def use_accession(data):
        source = data["regions"]["hc"][0]
        source.pop("path")
        source["accession"] = "NM_123456.2"

    update_manifest(manifest, use_accession)
    result = build(manifest, fake_builder, tmp_path)
    url = opener.requests[0][0]
    assert urlsplit(url).netloc == "eutils.ncbi.nlm.nih.gov"
    query = parse_qs(urlsplit(url).query)
    assert query["db"] == ["nuccore"]
    assert query["id"] == ["NM_123456.2"]
    assert query["rettype"] == ["fasta"]
    assert query["retmode"] == ["text"]
    assert result["manifest"]["regions"]["hc"]["sources"][0]["accession"] == "NM_123456.2"


@pytest.mark.parametrize("change, message", [
    ({"orientation": "reverse"}, "orientation"),
    ({"orientation": []}, "orientation"),
    ({"orientation": {}}, "orientation"),
    ({"locus": "TRB"}, "Invalid locus"),
    ({"region": "C"}, "Source region"),
    ({"url": "https://example.org/a.fa"}, "exactly one"),
])
def test_invalid_source_declarations_are_rejected(manifest, fake_builder, tmp_path, change, message):
    update_manifest(manifest, lambda data: data["regions"]["hv"][0].update(change))
    with pytest.raises(builder.ReferenceBuildError, match=message):
        build(manifest, fake_builder, tmp_path)


def test_variable_and_constant_locus_sets_must_match(manifest, fake_builder, tmp_path):
    update_manifest(manifest, lambda data: data["regions"]["lc"][0].update(locus="IGL"))
    with pytest.raises(builder.ReferenceBuildError, match="loci must match"):
        build(manifest, fake_builder, tmp_path)


@pytest.mark.parametrize("fasta,message", [
    (">x\nACGT\n>x\nTGCA\n", "Duplicate FASTA"),
    (">x\n", "Empty FASTA"),
    ("ACGT\n", "precedes a header"),
    (">x\nACGTZ\n", "Invalid FASTA sequence"),
    (">IGLV1\nACGT\n", "conflicts with declared locus"),
])
def test_invalid_fasta_is_rejected(manifest, fake_builder, tmp_path, fasta, message):
    (manifest.parent / "hv.fa").write_text(fasta)
    with pytest.raises(builder.ReferenceBuildError, match=message):
        build(manifest, fake_builder, tmp_path)


def test_tcr_mixed_loci_can_use_nonhuman_channel_layout(manifest, fake_builder, tmp_path):
    def tcr(data):
        for region, sources in data["regions"].items():
            sources[0]["locus"] = "TRB" if region.startswith("h") else "TRA"
        for region in ("hv", "hc"):
            extra = dict(data["regions"][region][0], locus="TRD")
            data["regions"][region].append(extra)

    update_manifest(manifest, tcr)
    result = builder.build_references(manifest, tmp_path / "references", "example_species", "TCR", fake_builder)
    assert result["manifest"]["regions"]["hv"]["loci"] == ["TRB", "TRD"]
    assert result["manifest"]["regions"]["lc"]["loci"] == ["TRA"]


def test_download_size_and_insecure_redirect_are_rejected(monkeypatch):
    monkeypatch.setattr(builder, "build_opener", lambda *_: FakeOpener(b"x" * 11))
    with pytest.raises(builder.ReferenceBuildError, match="max_bytes"):
        builder._Downloader(30, 10).fetch("https://example.org/file.fa")
    with pytest.raises(builder.ReferenceBuildError, match="HTTPS"):
        builder._HTTPSRedirectHandler().redirect_request(
            Request("https://example.org/file.fa"), None, 302, "redirect", {}, "http://example.org/file.fa")


def test_truncated_download_is_rejected_even_if_prefix_is_valid_fasta(monkeypatch):
    class Truncated:
        def open(self, request, **_kwargs):
            return FakeResponse(b">valid\nACGT\n", request.full_url, {"Content-Length": "100"})

    monkeypatch.setattr(builder, "build_opener", lambda *_: Truncated())
    with pytest.raises(builder.ReferenceBuildError, match="Truncated reference download"):
        builder._Downloader(30, 1000).fetch("https://example.org/truncated.fa")


def test_species_path_escape_is_rejected(manifest, fake_builder, tmp_path):
    with pytest.raises(builder.ReferenceBuildError, match="lowercase slug"):
        builder.build_references(manifest, tmp_path, "../outside", "BCR", fake_builder)


def test_cli_build_emits_json(manifest, fake_builder, tmp_path, capsys):
    assert builder.main(["build", "--manifest", str(manifest), "--output", str(tmp_path / "references"),
                         "--species", "example_species", "--receptor", "BCR",
                         "--bowtie-build", str(fake_builder)]) == 0
    assert json.loads(capsys.readouterr().out)["manifest"]["species"] == "example_species"


class SearchOpener:
    def __init__(self, identifiers=None, total=42):
        self.identifiers = ["11", "22"] if identifiers is None else identifiers
        self.total = total
        self.requests = []

    def open(self, request, timeout):
        self.requests.append(request.full_url)
        if "esearch.fcgi" in request.full_url:
            data = {"esearchresult": {"count": str(self.total), "idlist": self.identifiers,
                                       "querytranslation": "translated organism query"}}
        else:
            data = {"result": {"uids": self.identifiers,
                               **{uid: {"caption": "NM_" + uid, "accessionversion": "NM_" + uid + ".2",
                                        "title": "Candidate immunoglobulin transcript " + uid,
                                        "slen": 900 + int(uid), "taxid": 7955}
                                  for uid in self.identifiers}}}
        return FakeResponse(json.dumps(data).encode(), request.full_url)


def test_organism_search_returns_bounded_unclassified_candidates(monkeypatch):
    opener = SearchOpener()
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)
    monkeypatch.setattr(builder.time, "sleep", lambda _: None)
    result = builder.search_references("Danio rerio", "immunoglobulin", limit=2)
    assert result["total_matches"] == 42
    assert result["returned_count"] == 2
    assert result["candidates"][0]["accession_version"] == "NM_11.2"
    assert result["candidates"][0]["length"] == 911
    assert result["candidates"][0]["taxid"] == 7955
    assert not result["automatic_locus_assignment"]
    assert "locus" not in result["candidates"][0]
    search_query = parse_qs(urlsplit(opener.requests[0]).query)
    assert search_query["term"] == ['"Danio rerio"[Organism] AND (immunoglobulin)']
    assert search_query["retmax"] == ["2"]
    summary_query = parse_qs(urlsplit(opener.requests[1]).query)
    assert summary_query["id"] == ["11,22"]
    assert summary_query["version"] == ["2.0"]
    assert len(result["response_sha256"]) == 2


def test_empty_search_does_not_fetch_summaries(monkeypatch):
    opener = SearchOpener([], total=0)
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)
    result = builder.search_references("Danio rerio", "missing transcript")
    assert result["candidates"] == []
    assert len(opener.requests) == 1


def test_search_rejects_oversized_server_response(monkeypatch):
    opener = SearchOpener()
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)
    with pytest.raises(builder.ReferenceBuildError, match="oversized identifier"):
        builder.search_references("Danio rerio", "immunoglobulin", limit=1)


def test_http_failure_is_clear(monkeypatch):
    class Unavailable:
        def open(self, *_args, **_kwargs):
            raise URLError("service unavailable")

    monkeypatch.setattr(builder, "build_opener", lambda *_: Unavailable())
    with pytest.raises(builder.ReferenceBuildError, match="Cannot download.*service unavailable"):
        builder.search_references("Danio rerio", "immunoglobulin")


@pytest.mark.parametrize("limit", [0, 101, True])
def test_search_limit_is_checked_before_network(limit):
    with pytest.raises(builder.ReferenceBuildError, match="Search limit"):
        builder.search_references("Danio rerio", "immunoglobulin", limit=limit)


def test_search_cli_writes_candidate_json(tmp_path, monkeypatch):
    opener = SearchOpener()
    monkeypatch.setattr(builder, "build_opener", lambda *_: opener)
    monkeypatch.setattr(builder.time, "sleep", lambda _: None)
    output = tmp_path / "search" / "candidates.json"
    assert builder.main(["search", "--organism", "Danio rerio", "--query", "immunoglobulin",
                         "--limit", "2", "--output", str(output)]) == 0
    assert json.loads(output.read_text())["returned_count"] == 2


@pytest.mark.parametrize("update,message", [
    (lambda data: data.update(unknown=True), "top-level keys"),
    (lambda data: data.update(schema_version=2), "schema_version"),
    (lambda data: data["regions"].pop("lc"), "exactly hv"),
    (lambda data: data["regions"].update(hv=[]), "nonempty list"),
    (lambda data: data["regions"]["hv"][0].update(guess_locus=True), "unknown source field"),
    (lambda data: data["regions"]["hv"][0].update(path=""), "nonempty string"),
    (lambda data: data["regions"]["hv"][0].update(path="invalid\x00path.fa"), "NUL character"),
    (lambda data: data["regions"]["hv"][0].update(sha256="not a checksum"), "64 hexadecimal"),
    (lambda data: data.update(download={"retries": 4}), "only timeout and max_bytes"),
    (lambda data: data.update(download={"timeout": 0}), "Download timeout"),
    (lambda data: data.update(download={"max_bytes": 0}), "max_bytes"),
    (lambda data: data.update(metadata={"unsupported": "value"}), "metadata may contain"),
    (lambda data: data.update(metadata={"scientific_name": 123}), "nonempty string"),
    (lambda data: data.update(metadata={"taxon_id": True}), "positive integer"),
])
def test_manifest_schema_rejects_unusable_or_ambiguous_declarations(manifest, update, message):
    update_manifest(manifest, update)
    with pytest.raises(builder.ReferenceBuildError, match=message):
        builder._read_manifest(manifest, "BCR")


@pytest.mark.parametrize("raw,message", [
    (b'{"regions":{},"regions":{}}', "Duplicate JSON key"),
    (b'{not valid json}', "Cannot read reference manifest"),
    (b" " * 2_000_001, "exceeds 2 MB"),
])
def test_manifest_file_errors_are_explicit(tmp_path, raw, message):
    path = tmp_path / "manifest.json"
    path.write_bytes(raw)
    with pytest.raises(builder.ReferenceBuildError, match=message):
        builder._read_manifest(path, "BCR")


def test_missing_and_oversized_local_sources_do_not_install(manifest, fake_builder, tmp_path):
    (manifest.parent / "hv.fa").unlink()
    with pytest.raises(builder.ReferenceBuildError, match="not a regular file"):
        build(manifest, fake_builder, tmp_path)
    (manifest.parent / "hv.fa").write_text(">gene\nACGTACGTACGT\n")
    update_manifest(manifest, lambda data: data.update(download={"max_bytes": 10}))
    with pytest.raises(builder.ReferenceBuildError, match="exceeds max_bytes"):
        build(manifest, fake_builder, tmp_path)
    assert not (tmp_path / "references/example_species/BCR").exists()


def test_duplicate_identifiers_across_sources_are_rejected(manifest, fake_builder, tmp_path):
    update_manifest(manifest, lambda data: data["regions"]["hv"].append(dict(data["regions"]["hv"][0])))
    with pytest.raises(builder.ReferenceBuildError, match="Duplicate merged FASTA identifier"):
        build(manifest, fake_builder, tmp_path)


@pytest.mark.parametrize("content,message", [
    (b"\xff", "uncompressed ASCII FASTA"),
    (b">   \nACGT\n", "Invalid FASTA identifier"),
    (b"\n\n", "no FASTA records"),
])
def test_empty_or_non_fasta_payload_is_rejected(content, message):
    with pytest.raises(builder.ReferenceBuildError, match=message):
        builder._parse_fasta(content, {"locus": "IGH", "orientation": "forward"}, "hv")


def test_existing_build_lock_and_symlink_targets_are_preserved(manifest, fake_builder, tmp_path):
    parent = tmp_path / "references/example_species"
    parent.mkdir(parents=True)
    lock = parent / ".BCR.reference-build.lock"
    lock.write_text("other build")
    with pytest.raises(builder.ReferenceBuildError, match="Another build owns"):
        build(manifest, fake_builder, tmp_path)
    assert lock.read_text() == "other build"
    lock.unlink()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep").write_text("untouched")
    (parent / "BCR").symlink_to(outside, target_is_directory=True)
    with pytest.raises(builder.ReferenceBuildError, match="not a file or symlink"):
        build(manifest, fake_builder, tmp_path, overwrite=True)
    assert (outside / "keep").read_text() == "untouched"


def test_builder_directory_is_supported_and_missing_tool_is_clear(fake_builder, tmp_path):
    folder = tmp_path / "tools"
    folder.mkdir()
    (folder / "bowtie2-build").symlink_to(fake_builder)
    executable, version = builder._resolve_builder(folder)
    assert Path(executable).resolve() == fake_builder.resolve()
    assert version == "test bowtie2-build 1.0"
    with pytest.raises(builder.ReferenceBuildError, match="not executable"):
        builder._resolve_builder(tmp_path / "missing-builder")


@pytest.mark.parametrize("url", [
    "http://example.org/source.fa", "https://name:secret@example.org/a.fa",
    "https://example.org:invalid/a.fa", "https://example.org:99999/a.fa",
    "https://[broken/a.fa", 17,
])
def test_unsafe_or_malformed_source_urls_are_rejected(url):
    with pytest.raises(builder.ReferenceBuildError):
        builder._https_url(url)


def test_safe_https_redirect_is_followed():
    request = builder._HTTPSRedirectHandler().redirect_request(
        Request("https://example.org/old.fa"), None, 302, "redirect", {}, "https://example.org/new.fa")
    assert request.full_url == "https://example.org/new.fa"


@pytest.mark.parametrize("length,message", [("unknown", "Invalid Content-Length"),
                                             ("1001", "exceeds max_bytes"), ("-1", "exceeds max_bytes")])
def test_download_length_metadata_is_validated(monkeypatch, length, message):
    class ResponseOpener:
        def open(self, request, **kwargs):
            return FakeResponse(b">x\nACGT\n", request.full_url, {"Content-Length": length})

    monkeypatch.setattr(builder, "build_opener", lambda *_: ResponseOpener())
    with pytest.raises(builder.ReferenceBuildError, match=message):
        builder._Downloader(30, 1000).fetch("https://example.org/source.fa")


@pytest.mark.parametrize("delay_during", ["open", "read"])
def test_download_elapsed_deadline_is_enforced(monkeypatch, delay_during):
    clock = [0]

    class SlowResponse(FakeResponse):
        def read1(self, size):
            if delay_during == "read":
                clock[0] = 31
            return super().read1(size)

    class SlowOpener:
        def open(self, request, **kwargs):
            if delay_during == "open":
                clock[0] = 31
            return SlowResponse(b">x\nACGT\n", request.full_url)

    monkeypatch.setattr(builder.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(builder, "build_opener", lambda *_: SlowOpener())
    with pytest.raises(builder.ReferenceBuildError, match="exceeded its timeout"):
        builder._Downloader(30, 1000).fetch("https://example.org/source.fa")


@pytest.mark.parametrize("payload,message", [
    (b"<html>temporary outage</html>", "invalid JSON"),
    (b'{"error":"API rate limit exceeded"}', "esearch.fcgi failed"),
    (b'{"esearchresult":{"errorlist":{"phrasesnotfound":["bad taxon"]}}}', "nucleotide search failed"),
    (b'{"esearchresult":{"count":"invalid","idlist":[]}}', "invalid result count"),
    (b'{"esearchresult":{"count":"0","idlist":["11"]}}', "count is inconsistent"),
])
def test_ncbi_protocol_errors_do_not_become_empty_success(monkeypatch, payload, message):
    monkeypatch.setattr(builder, "build_opener", lambda *_: FakeOpener(payload))
    with pytest.raises(builder.ReferenceBuildError, match=message):
        builder.search_references("Danio rerio", "immunoglobulin")


def test_search_cli_stdout_and_error_preserve_existing_output(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(builder, "build_opener", lambda *_: SearchOpener([], total=0))
    assert builder.main(["search", "--organism", "Danio rerio", "--query", "immunoglobulin"]) == 0
    assert json.loads(capsys.readouterr().out)["returned_count"] == 0
    output = tmp_path / "previous.json"
    output.write_text("previous result")
    with pytest.raises(SystemExit) as error:
        builder.main(["search", "--organism", "Danio rerio", "--query", "immunoglobulin",
                      "--limit", "0", "--output", str(output)])
    assert error.value.code == 2
    assert "basic-references: Search limit" in capsys.readouterr().err
    assert output.read_text() == "previous result"
