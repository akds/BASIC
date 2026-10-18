"""Retrieve declared receptor references and install checked Bowtie2 indexes.

Retrieval is explicit by local file, HTTPS URL, or NCBI nucleotide accession.
This module does not infer germline annotations from a species name.
"""

from __future__ import annotations

import argparse
import hashlib
import hmac
from http.client import HTTPException
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener
import uuid

from . import __version__
from .alignment import AlignmentError, validate_indexes
from .assembly import normalize_sequence, revcom


REGIONS = ("hv", "hc", "lv", "lc")
LOCI = frozenset(("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"))
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_MAX_BYTES = 50_000_000


class ReferenceBuildError(RuntimeError):
    """A declared source or reference build could not be completed safely."""


def _object_without_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ReferenceBuildError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _https_url(url):
    if not isinstance(url, str):
        raise ReferenceBuildError("Source URL must be a string")
    try:
        parsed = urlsplit(url)
        # Accessing .port validates both nonnumeric and out-of-range ports.
        _ = parsed.port
    except ValueError as error:
        raise ReferenceBuildError(f"Invalid reference URL: {error}") from error
    if (parsed.scheme != "https" or not parsed.hostname or parsed.username is not None
            or parsed.password is not None or parsed.fragment
            or any(character.isspace() for character in url)):
        raise ReferenceBuildError("Sources and redirects must use HTTPS without credentials or fragments")
    return url


class _HTTPSRedirectHandler(HTTPRedirectHandler):
    max_redirections = 5

    def redirect_request(self, request, response, code, message, headers, newurl):
        _https_url(newurl)
        return super().redirect_request(request, response, code, message, headers, newurl)


def accession_url(accession):
    """Use NCBI's documented nucleotide EFetch FASTA interface."""
    if not isinstance(accession, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*(?:\.[0-9]+)?", accession):
        raise ReferenceBuildError("accession must be one nucleotide accession, optionally with a version")
    return NCBI_EFETCH + "?" + urlencode({
        "db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
        "tool": "BASIC-reference-builder"})


class _Downloader:
    def __init__(self, timeout, max_bytes):
        self.timeout = timeout
        self.max_bytes = max_bytes
        self.opener = build_opener(_HTTPSRedirectHandler())
        self.last_ncbi_request = None

    def fetch(self, url):
        _https_url(url)
        if urlsplit(url).hostname == "eutils.ncbi.nlm.nih.gov":
            if self.last_ncbi_request is not None:
                time.sleep(max(0.0, 0.35 - (time.monotonic() - self.last_ncbi_request)))
            self.last_ncbi_request = time.monotonic()
        request = Request(url, headers={"User-Agent": f"BASIC/{__version__} reference-builder",
                                        "Accept": "text/plain", "Accept-Encoding": "identity"})
        started = time.monotonic()
        try:
            with self.opener.open(request, timeout=self.timeout) as response:
                final_url = _https_url(response.geturl())
                declared_length = response.headers.get("Content-Length")
                if declared_length is not None:
                    try:
                        length = int(declared_length)
                    except ValueError as error:
                        raise ReferenceBuildError("Invalid Content-Length in reference download") from error
                    if length < 0 or length > self.max_bytes:
                        raise ReferenceBuildError("Reference download exceeds max_bytes")
                chunks, received = [], 0
                read = getattr(response, "read1", response.read)
                while True:
                    if time.monotonic() - started > self.timeout:
                        raise ReferenceBuildError("Reference download exceeded its timeout")
                    chunk = read(min(65536, self.max_bytes + 1 - received))
                    if time.monotonic() - started > self.timeout:
                        raise ReferenceBuildError("Reference download exceeded its timeout")
                    if not chunk:
                        break
                    chunks.append(chunk)
                    received += len(chunk)
                    if received > self.max_bytes:
                        raise ReferenceBuildError("Reference download exceeds max_bytes")
                if declared_length is not None and received != length:
                    raise ReferenceBuildError(
                        f"Truncated reference download: Content-Length was {length}, received {received}")
                return b"".join(chunks), final_url
        except (HTTPError, URLError, HTTPException, TimeoutError, OSError) as error:
            raise ReferenceBuildError(f"Cannot download {url}: {error}") from error


def _read_manifest(path, receptor):
    try:
        with path.open("rb") as handle:
            raw = handle.read(2_000_001)
        if len(raw) > 2_000_000:
            raise ReferenceBuildError("Reference manifest exceeds 2 MB")
        data = json.loads(raw, object_pairs_hook=_object_without_duplicates)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ReferenceBuildError(f"Cannot read reference manifest {path}: {error}") from error
    if not isinstance(data, dict) or set(data) - {"schema_version", "regions", "download", "metadata"}:
        raise ReferenceBuildError("Manifest must contain regions and only documented top-level keys")
    if type(data.get("schema_version", 1)) is not int or data.get("schema_version", 1) != 1:
        raise ReferenceBuildError("Only reference manifest schema_version 1 is supported")
    regions = data.get("regions")
    if not isinstance(regions, dict) or set(regions) != set(REGIONS):
        raise ReferenceBuildError("Manifest regions must contain exactly hv, hc, lv, and lc")
    for region, sources in regions.items():
        if not isinstance(sources, list) or not sources:
            raise ReferenceBuildError(f"{region} must contain a nonempty list of sources")
        for source in sources:
            if not isinstance(source, dict) or set(source) - {
                    "path", "url", "accession", "orientation", "locus", "region", "sha256"}:
                raise ReferenceBuildError(f"Invalid or unknown source field in {region}")
            kinds = set(source) & {"path", "url", "accession"}
            if len(kinds) != 1:
                raise ReferenceBuildError(f"Every {region} source needs exactly one path, url, or accession")
            kind = next(iter(kinds))
            if not isinstance(source[kind], str) or not source[kind].strip():
                raise ReferenceBuildError(f"Source {kind} must be a nonempty string")
            if "\x00" in source[kind]:
                raise ReferenceBuildError(f"Source {kind} contains a NUL character")
            if not isinstance(source.get("orientation"), str) or source["orientation"] not in {"forward", "as_is"}:
                raise ReferenceBuildError(f"Every {region} source must declare orientation forward or as_is")
            expected_region = "V" if region.endswith("v") else "C"
            if source.get("region", expected_region) != expected_region:
                raise ReferenceBuildError(f"Source region must be {expected_region} inside {region}")
            locus = source.get("locus")
            allowed = ({"IGH"} if region.startswith("h") else {"IGK", "IGL"}) if receptor == "BCR" else {"TRA", "TRB", "TRG", "TRD"}
            if not isinstance(locus, str) or locus not in allowed:
                raise ReferenceBuildError(f"Invalid locus for {receptor} {region}; choose from {', '.join(sorted(allowed))}")
            if "sha256" in source and (not isinstance(source["sha256"], str)
                    or not re.fullmatch(r"[0-9a-fA-F]{64}", source["sha256"])):
                raise ReferenceBuildError("sha256 must contain exactly 64 hexadecimal characters")
            if kind == "url":
                _https_url(source[kind])
            if kind == "accession":
                accession_url(source[kind])
    for variable, constant in (("hv", "hc"), ("lv", "lc")):
        if {item["locus"] for item in regions[variable]} != {item["locus"] for item in regions[constant]}:
            raise ReferenceBuildError(f"Declared loci must match between {variable} and {constant}")
    download = data.get("download", {})
    if not isinstance(download, dict) or set(download) - {"timeout", "max_bytes"}:
        raise ReferenceBuildError("download may contain only timeout and max_bytes")
    timeout = download.get("timeout", 30)
    maximum = download.get("max_bytes", DEFAULT_MAX_BYTES)
    if (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
            or not math.isfinite(timeout) or not 0 < timeout <= 300):
        raise ReferenceBuildError("Download timeout must be greater than 0 and at most 300 seconds")
    if isinstance(maximum, bool) or not isinstance(maximum, int) or not 1 <= maximum <= 1_000_000_000:
        raise ReferenceBuildError("max_bytes must be an integer from 1 to 1000000000")
    metadata = data.get("metadata", {})
    if not isinstance(metadata, dict) or set(metadata) - {"scientific_name", "taxon_id", "notes"}:
        raise ReferenceBuildError("metadata may contain only scientific_name, taxon_id, and notes")
    for key in ("scientific_name", "notes"):
        if key in metadata and (not isinstance(metadata[key], str) or not metadata[key].strip()):
            raise ReferenceBuildError(f"metadata.{key} must be a nonempty string")
    if "taxon_id" in metadata and (isinstance(metadata["taxon_id"], bool)
            or not isinstance(metadata["taxon_id"], int) or metadata["taxon_id"] <= 0):
        raise ReferenceBuildError("metadata.taxon_id must be a positive integer")
    return data, raw, float(timeout), maximum


def _source_data(source, manifest_directory, downloader):
    if "path" in source:
        path = Path(source["path"]).expanduser()
        if not path.is_absolute():
            path = manifest_directory / path
        path = path.resolve()
        try:
            if not path.is_file():
                raise ReferenceBuildError(f"Reference source is not a regular file: {path}")
            with path.open("rb") as handle:
                content = handle.read(downloader.max_bytes + 1)
        except OSError as error:
            raise ReferenceBuildError(f"Cannot read reference source {path}: {error}") from error
        if len(content) > downloader.max_bytes:
            raise ReferenceBuildError(f"Reference source exceeds max_bytes: {path}")
        provenance = {"type": "path", "path": str(path), "declared_path": source["path"]}
    else:
        url = source.get("url") or accession_url(source["accession"])
        content, final_url = downloader.fetch(url)
        provenance = {"type": "url" if "url" in source else "accession", "url": url,
                      "resolved_url": final_url}
        if "accession" in source:
            provenance["accession"] = source["accession"]
    digest = hashlib.sha256(content).hexdigest()
    if "sha256" in source and not hmac.compare_digest(digest, source["sha256"].lower()):
        raise ReferenceBuildError(f"SHA256 mismatch for source {provenance}")
    provenance["sha256"] = digest
    if "sha256" in source:
        provenance["expected_sha256"] = source["sha256"].lower()
    return content, provenance


def _parse_fasta(content, source, region):
    try:
        text = content.decode("ascii")
    except UnicodeDecodeError as error:
        raise ReferenceBuildError("References must be uncompressed ASCII FASTA") from error
    records, seen = [], set()
    identifier, description, parts = None, "", []

    def save():
        if identifier is None:
            return
        try:
            sequence = normalize_sequence("".join(parts))
        except ValueError as error:
            raise ReferenceBuildError(f"Invalid FASTA sequence {identifier}: {error}") from error
        if not sequence:
            raise ReferenceBuildError(f"Empty FASTA sequence: {identifier}")
        if identifier in seen:
            raise ReferenceBuildError(f"Duplicate FASTA identifier: {identifier}")
        seen.add(identifier)
        named_loci = set(re.findall(r"(?:^|[|:_])(IGH|IGK|IGL|TRA|TRB|TRG|TRD)", identifier.upper()))
        if named_loci - {source["locus"]}:
            raise ReferenceBuildError(f"FASTA identifier {identifier} conflicts with declared locus {source['locus']}")
        if region.endswith("c") and source["orientation"] == "forward":
            sequence = revcom(sequence)
        prefix = source["locus"] + "|"
        oriented_id = identifier if identifier.startswith(prefix) else prefix + identifier
        records.append((oriented_id, description, sequence))

    for number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            save()
            words = line[1:].split(maxsplit=1)
            if not words or any(ord(char) < 33 or ord(char) > 126 for char in words[0]):
                raise ReferenceBuildError(f"Invalid FASTA identifier on line {number}")
            identifier = words[0]
            description = words[1] if len(words) == 2 else ""
            parts = []
        else:
            if identifier is None:
                raise ReferenceBuildError(f"FASTA sequence precedes a header on line {number}")
            parts.append(line)
    save()
    if not records:
        raise ReferenceBuildError("Reference source contains no FASTA records")
    return records


def _run_builder(command, timeout):
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   start_new_session=(os.name == "posix"))
    except OSError as error:
        raise ReferenceBuildError(f"Cannot start Bowtie2 builder: {error}") from error
    try:
        output, _ = process.communicate(timeout=timeout)
    except BaseException:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.communicate()
        raise
    return process.returncode, output


def _resolve_builder(requested):
    path = Path(requested).expanduser()
    if path.is_dir():
        path /= "bowtie2-build"
    resolved = str(path.resolve()) if path.is_file() else shutil.which(str(path))
    if not resolved or not os.access(resolved, os.X_OK):
        raise ReferenceBuildError(f"Bowtie2 builder is not executable: {requested}")
    code, output = _run_builder([resolved, "--version"], 30)
    lines = output.decode("utf-8", errors="replace").strip().splitlines()
    if code or not lines:
        raise ReferenceBuildError("bowtie2-build --version failed or returned no version")
    return resolved, lines[0]


def _install(staging, target, overwrite):
    if target.is_symlink() or (target.exists() and not target.is_dir()):
        raise ReferenceBuildError(f"Target must be a directory, not a file or symlink: {target}")
    if not target.exists():
        os.replace(staging, target)
        return None
    if not overwrite:
        raise ReferenceBuildError(f"Reference target already exists; use --overwrite: {target}")
    backup = target.with_name(f".{target.name}.backup-{uuid.uuid4().hex}")
    os.replace(target, backup)
    try:
        os.replace(staging, target)
    except BaseException as install_error:
        try:
            os.replace(backup, target)
        except OSError as rollback_error:
            raise ReferenceBuildError(
                f"Installation and rollback failed; previous references remain at {backup}: {rollback_error}") from install_error
        raise
    try:
        shutil.rmtree(backup)
    except OSError:
        return str(backup)
    return None


def build_references(manifest_path, output_root, species, receptor,
                     bowtie_build="bowtie2-build", overwrite=False, *, build_timeout=600):
    """Stage all four indexes, validate them, then install with rollback.

    The target is ``output_root/species/receptor``. Overwrite swaps directories
    with a rollback backup; readers should not run during that brief rename
    window. Source/download/build failures never modify an existing target.
    """
    if not isinstance(species, str) or not re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", species):
        raise ReferenceBuildError("species must be a lowercase slug using letters, digits, _ or -")
    if not isinstance(receptor, str) or receptor not in {"BCR", "TCR"}:
        raise ReferenceBuildError("receptor must be BCR or TCR")
    if (isinstance(build_timeout, bool) or not isinstance(build_timeout, (int, float))
            or not math.isfinite(build_timeout) or build_timeout <= 0):
        raise ReferenceBuildError("build_timeout must be finite and positive")
    manifest_path = Path(manifest_path).expanduser().resolve()
    data, raw_manifest, timeout, maximum = _read_manifest(manifest_path, receptor)
    root = Path(output_root).expanduser().resolve()
    parent = root / species
    if parent.is_symlink():
        raise ReferenceBuildError("Species destination must not be a symlink")
    parent.mkdir(parents=True, exist_ok=True)
    target = parent / receptor
    lock = parent / f".{receptor}.reference-build.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ReferenceBuildError(f"Another build owns {lock}; inspect the owner before removing a stale lock") from error
    staging = None
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(str(os.getpid()) + "\n")
        if target.is_symlink() or (target.exists() and not target.is_dir()):
            raise ReferenceBuildError(f"Target must be a directory, not a file or symlink: {target}")
        if target.exists() and not overwrite:
            raise ReferenceBuildError(f"Reference target already exists; use --overwrite: {target}")
        executable, version = _resolve_builder(bowtie_build)
        staging = Path(tempfile.mkdtemp(prefix=f".{receptor}.building-", dir=parent))
        downloader = _Downloader(timeout, maximum)
        provenance = {
            "schema_version": 1, "species": species, "receptor": receptor,
            "metadata": data.get("metadata", {}), "basic_version": __version__,
            "source_manifest": {"path": str(manifest_path),
                                "sha256": hashlib.sha256(raw_manifest).hexdigest()},
            "builder": {"executable": executable, "version": version,
                        "timeout_seconds": build_timeout},
            "download_limits": {"timeout_seconds": timeout, "max_bytes": maximum},
            "regions": {},
        }
        for region in REGIONS:
            identifiers, records, sources = set(), [], []
            for source in data["regions"][region]:
                content, source_provenance = _source_data(source, manifest_path.parent, downloader)
                fetched = _parse_fasta(content, source, region)
                for record in fetched:
                    if record[0] in identifiers:
                        raise ReferenceBuildError(f"Duplicate merged FASTA identifier in {region}: {record[0]}")
                    identifiers.add(record[0])
                    records.append(record)
                source_provenance.update({"orientation": source["orientation"], "locus": source["locus"],
                                          "region": "V" if region.endswith("v") else "C",
                                          "reverse_complemented": region.endswith("c") and source["orientation"] == "forward",
                                          "records": len(fetched)})
                sources.append(source_provenance)
            fasta = "".join(f">{identifier}{' ' + description if description else ''}\n{sequence}\n"
                            for identifier, description, sequence in records).encode("ascii")
            (staging / region).write_bytes(fasta)
            provenance["regions"][region] = {"fasta": region, "sha256": hashlib.sha256(fasta).hexdigest(),
                                             "records": len(records), "loci": sorted({s["locus"] for s in sources}),
                                             "sources": sources}
        for region in REGIONS:
            code, output = _run_builder([executable, str(staging / region), str(staging / region)], build_timeout)
            (staging / f"build_{region}.log").write_bytes(output)
            if code:
                detail = output[-2000:].decode("utf-8", errors="replace").strip()
                raise ReferenceBuildError(f"bowtie2-build failed for {region} (exit {code}): {detail}")
        validate_indexes(staging)
        (staging / "reference_manifest.json").write_text(
            json.dumps(provenance, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        backup = _install(staging, target, overwrite)
        return {"path": str(target), "manifest": provenance, "retained_backup": backup}
    except (OSError, AlignmentError, subprocess.TimeoutExpired) as error:
        raise ReferenceBuildError(f"Reference build failed: {error}") from error
    finally:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)
        lock.unlink(missing_ok=True)


def search_references(organism, query, limit=10, *, timeout=30, max_bytes=5_000_000):
    """Find bounded, uncurated NCBI nucleotide candidates for manual selection.

    Search never assigns V/C regions or loci and never installs a reference.
    A candidate title is database metadata, not evidence of germline status.
    """
    if (not isinstance(organism, str) or not organism.strip() or len(organism) > 200
            or any(character in '\"[]' or ord(character) < 32 for character in organism)):
        raise ReferenceBuildError("organism must be a nonempty name without query delimiters")
    if (not isinstance(query, str) or not query.strip() or len(query) > 1000
            or any(ord(character) < 32 for character in query)):
        raise ReferenceBuildError("query must be nonempty and at most 1000 characters")
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
        raise ReferenceBuildError("Search limit must be an integer from 1 to 100")
    if (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
            or not math.isfinite(timeout) or not 0 < timeout <= 300):
        raise ReferenceBuildError("Search timeout must be greater than 0 and at most 300 seconds")
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or not 1 <= max_bytes <= 1_000_000_000:
        raise ReferenceBuildError("Search max_bytes must be an integer from 1 to 1000000000")
    organism, query = organism.strip(), query.strip()
    term = f'"{organism}"[Organism] AND ({query})'
    downloader = _Downloader(timeout, max_bytes)
    urls, response_hashes = [], []

    def retrieve(endpoint, parameters):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint + "?" + urlencode({
            "db": "nuccore", "retmode": "json", "tool": "BASIC-reference-builder", **parameters})
        raw, final_url = downloader.fetch(url)
        urls.append(final_url)
        response_hashes.append(hashlib.sha256(raw).hexdigest())
        try:
            parsed = json.loads(raw)
        except (UnicodeError, json.JSONDecodeError) as error:
            raise ReferenceBuildError(f"NCBI {endpoint} returned invalid JSON") from error
        if not isinstance(parsed, dict) or parsed.get("error"):
            raise ReferenceBuildError(f"NCBI {endpoint} failed: {parsed}")
        return parsed

    search = retrieve("esearch.fcgi", {"term": term, "retmax": limit}).get("esearchresult")
    if not isinstance(search, dict) or search.get("error") or search.get("errorlist"):
        raise ReferenceBuildError(f"NCBI nucleotide search failed: {search}")
    identifiers = search.get("idlist")
    if (not isinstance(identifiers, list) or len(identifiers) > limit
            or any(not isinstance(uid, str) or not uid.isdecimal() for uid in identifiers)
            or len(identifiers) != len(set(identifiers))):
        raise ReferenceBuildError("NCBI search returned an invalid or oversized identifier list")
    try:
        total = int(search["count"])
    except (KeyError, TypeError, ValueError) as error:
        raise ReferenceBuildError("NCBI search returned an invalid result count") from error
    if total < len(identifiers):
        raise ReferenceBuildError("NCBI search result count is inconsistent")
    candidates = []
    if identifiers:
        summaries = retrieve("esummary.fcgi", {"id": ",".join(identifiers), "version": "2.0"}).get("result")
        if not isinstance(summaries, dict):
            raise ReferenceBuildError("NCBI returned no nucleotide summaries")
        for uid in identifiers:
            record = summaries.get(uid)
            if not isinstance(record, dict) or record.get("error"):
                raise ReferenceBuildError(f"NCBI returned no usable summary for {uid}")
            accession_version = record.get("accessionversion") or None
            accession = record.get("caption") or accession_version
            title = record.get("title")
            if (not isinstance(accession, str) or not accession
                    or (accession_version is not None and not isinstance(accession_version, str))
                    or not isinstance(title, str)):
                raise ReferenceBuildError(f"NCBI summary is missing accession/title fields for {uid}")
            try:
                length = int(record["slen"])
                taxid = int(record.get("taxid", 0)) or None
            except (KeyError, TypeError, ValueError) as error:
                raise ReferenceBuildError(f"NCBI summary has invalid length/taxid for {uid}") from error
            if length < 1 or (taxid is not None and taxid < 1):
                raise ReferenceBuildError(f"NCBI summary has nonpositive length/taxid for {uid}")
            candidates.append({"uid": uid, "accession": accession,
                               "accession_version": accession_version, "title": title,
                               "length": length, "taxid": taxid,
                               "record_url": "https://www.ncbi.nlm.nih.gov/nuccore/" + uid})
    return {"schema_version": 1, "database": "nuccore", "organism": organism,
            "query": query, "search_term": term, "query_translation": search.get("querytranslation"),
            "limit": limit, "total_matches": total, "returned_count": len(candidates),
            "candidates": candidates, "candidate_status": "unvalidated_reference_candidates",
            "automatic_locus_assignment": False, "warnings": search.get("warninglist", {}),
            "source_urls": urls, "response_sha256": response_hashes,
            "download_limits": {"timeout_seconds": timeout, "max_bytes": max_bytes}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="Download/import declared FASTA and build all four indexes")
    build.add_argument("--manifest", required=True)
    build.add_argument("--output", required=True)
    build.add_argument("--species", required=True)
    build.add_argument("--receptor", required=True, choices=("BCR", "TCR"))
    build.add_argument("--bowtie-build", default="bowtie2-build")
    build.add_argument("--build-timeout", type=float, default=600)
    build.add_argument("--overwrite", action="store_true")
    search = commands.add_parser("search", help="Find NCBI nucleotide candidates by organism and query")
    search.add_argument("--organism", required=True)
    search.add_argument("--query", required=True)
    search.add_argument("--limit", type=int, default=10)
    search.add_argument("--timeout", type=float, default=30)
    search.add_argument("--max-bytes", type=int, default=5_000_000)
    search.add_argument("--output", type=Path, help="Write candidate JSON here instead of stdout")
    args = parser.parse_args(argv)
    try:
        if args.command == "search":
            result = search_references(args.organism, args.query, args.limit,
                                        timeout=args.timeout, max_bytes=args.max_bytes)
        else:
            result = build_references(args.manifest, args.output, args.species, args.receptor,
                                      args.bowtie_build, args.overwrite, build_timeout=args.build_timeout)
        serialized = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.command == "search" and args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=args.output.parent,
                                                 prefix=f".{args.output.name}.", delete=False) as handle:
                    temporary = Path(handle.name)
                    handle.write(serialized)
                os.replace(temporary, args.output)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        else:
            print(serialized, end="")
    except (ReferenceBuildError, OSError) as error:
        parser.exit(2, f"basic-references: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
