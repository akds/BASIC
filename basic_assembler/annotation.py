"""AIRR TSV export and optional, local IgBLAST annotation.

AIRR required columns are nullable; absence of an annotation is never evidence
for an unproductive receptor. IgBLAST provides biological calls without BASIC
guessing them from contig length or conserved motifs.
"""

import csv
import hashlib
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
import time

from .processes import cleanup_process


# AIRR Rearrangement schema 2.0 (also compatible with the 1.4 core).
REQUIRED_FIELDS = (
    "sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call",
    "j_call", "sequence_alignment", "germline_alignment", "junction",
    "junction_aa", "v_cigar", "d_cigar", "j_cigar",
)
ASSEMBLY_FIELDS = ("cell_id", "locus", "basic_status", "basic_ambiguous",
                   "basic_read_support")
BOOLEAN_FIELDS = frozenset(("rev_comp", "productive", "vj_in_frame", "stop_codon",
                            "complete_vdj", "basic_ambiguous"))
ALIASES = {"status": "basic_status", "ambiguous": "basic_ambiguous",
           "read_support": "basic_read_support"}


class AnnotationError(RuntimeError):
    """An annotation subprocess failed or produced unusable AIRR output."""


def _run_igblast(command, environment, timeout):
    """Reap the program and, on POSIX, its descendants when interrupted."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, encoding="utf-8", errors="replace", env=environment,
                               start_new_session=(os.name == "posix"))
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except BaseException:
        cleanup_process(process)
        raise
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def _cell(value, field):
    if value is None:
        return ""
    if isinstance(value, bool):
        value = "T" if value else "F"
    elif isinstance(value, (list, tuple, dict, set)):
        raise ValueError("AIRR field %s must be a scalar" % field)
    else:
        value = str(value)
    if any(char in value for char in "\t\r\n"):
        raise ValueError("AIRR field %s contains a tab or newline" % field)
    if field in BOOLEAN_FIELDS and value not in ("", "T", "F"):
        raise ValueError("AIRR boolean %s must be T, F, or null" % field)
    return value


def _normalize(record):
    return {ALIASES.get(key, key): _cell(value, ALIASES.get(key, key))
            for key, value in record.items()}


def _atomic_table(records, path, fields):
    destination = Path(path)
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=destination.parent,
                                         prefix=".%s." % destination.name,
                                         suffix=".tmp", delete=False) as handle:
            temp_path = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n", quoting=csv.QUOTE_NONE)
            writer.writeheader()
            writer.writerows(records)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, destination)
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()


def write_airr(records, path):
    """Atomically write records as AIRR TSV; unknown required values stay empty.

    ``status``, ``ambiguous``, and ``read_support`` map to ``basic_*`` columns.
    Support is deliberately not mapped to AIRR ``consensus_count``: BASIC's
    read-pool support need not equal the number of reads in a consensus.
    Additional scalar AIRR/custom columns are preserved. Returns the row count.
    """
    normalized = [_normalize(record) for record in records]
    fields = list(REQUIRED_FIELDS + ASSEMBLY_FIELDS)
    fields.extend(sorted(set().union(*(set(row) for row in normalized)) - set(fields)))
    for field in fields:
        _cell(field, "column name")
    identifiers = [row.get("sequence_id") for row in normalized if row.get("sequence_id")]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("AIRR sequence_id values must be unique")
    _atomic_table(normalized, path, fields)
    return len(normalized)


def read_airr(path):
    """Read AIRR TSV, validating the core header, row widths, and booleans.

    This is structural validation, not a complete AIRR schema validator.
    """
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        fields = reader.fieldnames or []
        missing = set(REQUIRED_FIELDS) - set(fields)
        if missing:
            raise AnnotationError("AIRR output is missing required columns: " + ", ".join(sorted(missing)))
        if len(fields) != len(set(fields)) or "" in fields:
            raise AnnotationError("AIRR output has duplicate or empty column names")
        records, seen = [], set()
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise AnnotationError("AIRR output has an invalid field count on line %d" % line)
            try:
                record = _normalize(row)
            except ValueError as error:
                raise AnnotationError(str(error)) from error
            identifier = record["sequence_id"]
            if identifier in seen or not identifier:
                raise AnnotationError("AIRR annotation requires unique, nonempty sequence_id values")
            seen.add(identifier)
            records.append(record)
    return records


def _fasta_ids(path):
    identifiers, active, has_sequence = [], None, False
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if active is not None and not has_sequence:
                    raise ValueError("Empty FASTA sequence: " + active)
                words = line[1:].split()
                if not words:
                    raise ValueError("FASTA identifier must not be empty")
                active, has_sequence = words[0], False
                identifiers.append(active)
            elif active is None:
                raise ValueError("FASTA sequence precedes its header")
            else:
                has_sequence = True
    if active is not None and not has_sequence:
        raise ValueError("Empty FASTA sequence: " + active)
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("FASTA identifiers must be unique for annotation")
    return identifiers


def _resolve_auxiliary(auxiliary_data, igdata):
    if not auxiliary_data:
        return None
    requested = Path(auxiliary_data).expanduser()
    candidates = [requested]
    data_root = igdata or os.environ.get("IGDATA")
    if data_root and not requested.is_absolute():
        candidates.append(Path(data_root).expanduser() / requested)
    for path in candidates:
        if path.is_file() and os.access(path, os.R_OK):
            return path.resolve()
    raise ValueError("IgBLAST auxiliary file is not readable: " + str(auxiliary_data))


def _explicit_prefix(value):
    path = Path(value).expanduser()
    return path.is_absolute() or "/" in str(value) or "\\" in str(value)


def _database_files(prefix):
    """Find nucleotide database files, including numbered volume shards."""
    prefix = Path(prefix)
    pattern = re.compile(re.escape(prefix.name) + r"(?:\.\d+)?\.n[a-z]{2}$")
    try:
        return sorted(path for path in prefix.parent.iterdir()
                      if pattern.fullmatch(path.name) and path.is_file())
    except (FileNotFoundError, NotADirectoryError):
        return []


def _find_database(value, relative_to=None):
    path = Path(value).expanduser()
    if path.is_absolute():
        candidates = [path]
    elif _explicit_prefix(value):
        candidates = [(relative_to or Path.cwd()) / path]
    else:
        roots = ([relative_to] if relative_to else []) + [Path.cwd()]
        roots += [Path(root).expanduser() for root in os.environ.get("BLASTDB", "").split(os.pathsep) if root]
        candidates = [root / path for root in roots]
    for candidate in candidates:
        files = _database_files(candidate)
        if files:
            return candidate.resolve(), files
    return None, []


def annotation_reference_fingerprints(*, germline_v, germline_j, germline_d=None,
                                      germline_c=None, auxiliary_data=None, igdata=None,
                                      organism="human", executable=None):
    """Fingerprint discoverable annotation references without inventing paths.

    Searches explicit prefixes, the working directory, and BLASTDB directories.
    Follows nucleotide alias DBLIST entries and common alias filter files.
    Hashes the selected organism's internal_data files and supplied auxiliary
    file. Does not read .ncbirc or claim to reproduce every BLAST lookup rule;
    missing resources, unsupported alias lists, cycles, and read failures appear
    in ``unresolved``. Discovered files are not a database-validity assertion.
    """
    found, unresolved, databases, visited = {}, [], [], set()

    def missing(role, requested, reason):
        unresolved.append({"role": role, "requested": str(requested), "reason": reason})

    def fingerprint(path, role):
        path = Path(path).resolve()
        key = str(path)
        if key in found:
            found[key]["roles"].add(role)
            return
        try:
            before = path.stat()
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                missing(role, path, "File changed while its checksum was computed")
                return
            found[key] = {"path": key, "bytes": after.st_size,
                          "sha256": digest.hexdigest(), "roles": {role}}
        except OSError as error:
            missing(role, path, "Cannot fingerprint file: " + str(error))

    def visit(value, role, relative_to=None, active=()):
        try:
            prefix, files = _find_database(value, relative_to)
        except OSError as error:
            missing(role, value, "Cannot inspect database directory: " + str(error))
            return
        if prefix is None:
            missing(role, value, "No nucleotide database files found in explicit path, working directory, or BLASTDB")
            return
        if prefix in active or len(active) >= 32:
            missing(role, prefix, "Cyclic alias or alias nesting exceeds 32 levels")
            return
        if (prefix, role) in visited:
            return
        visited.add((prefix, role))
        databases.append({"role": role, "requested": str(value), "resolved_prefix": str(prefix)})
        for path in files:
            fingerprint(path, role)
            if path.suffix != ".nal":
                continue
            try:
                # Alias files are small text records. A trailing backslash is
                # accepted as a continuation; shell commands are never run.
                text = path.read_text(encoding="utf-8").replace("\\\n", " ")
                for line in text.splitlines():
                    tokens = shlex.split(line, comments=True)
                    if not tokens:
                        continue
                    kind, values = tokens[0].upper(), tokens[1:]
                    if kind == "DBLIST":
                        for value in values:
                            visit(value, role, path.parent, active + (prefix,))
                    elif kind in {"GILIST", "OIDLIST", "TILIST", "SEQIDLIST", "TAXIDLIST"}:
                        for value in values:
                            dependency = Path(value).expanduser()
                            fingerprint(dependency if dependency.is_absolute() else path.parent / dependency, role)
                    elif kind.endswith("LIST"):
                        missing(role, path, "Alias list type not resolved: " + kind)
            except (OSError, UnicodeError, ValueError) as error:
                missing(role, path, "Cannot resolve alias file: " + str(error))

    for role, value in (("germline_v", germline_v), ("germline_d", germline_d),
                        ("germline_j", germline_j), ("germline_c", germline_c)):
        if value:
            visit(value, role)
    if auxiliary_data:
        try:
            fingerprint(_resolve_auxiliary(auxiliary_data, igdata), "auxiliary_data")
        except ValueError as error:
            missing("auxiliary_data", auxiliary_data, str(error))
    data_root = igdata or os.environ.get("IGDATA")
    if data_root:
        roots = [Path(data_root).expanduser()]
    else:
        roots = [Path.cwd()]
        if executable:
            executable_path = Path(executable).expanduser().resolve()
            roots.extend((executable_path.parent, executable_path.parent.parent))
    if Path(organism).name != organism or organism in (".", ".."):
        missing("internal_data", organism, "Organism must name a single internal_data subdirectory")
    else:
        directories = [root / "internal_data" / organism for root in roots]
        directory = next((path for path in directories if path.is_dir()), None)
        if directory is None:
            missing("internal_data", organism, "Selected organism's installation data could not be located")
        else:
            paths = [path for path in sorted(directory.rglob("*")) if path.is_file()]
            if not paths:
                missing("internal_data", directory, "Selected organism's installation directory contains no files")
            for path in paths:
                fingerprint(path, "internal_data")
    return {"files": [dict(item, roles=sorted(item["roles"])) for _, item in sorted(found.items())],
            "database_prefixes": databases, "unresolved": unresolved,
            "resolution_scope": "Explicit prefixes, cwd, BLASTDB, nucleotide aliases, supplied auxiliary file, and located organism internal_data; .ncbirc and implicit default databases are not resolved"}


def _validate_parameters(organism, receptor, germline_v, germline_j,
                         igdata, auxiliary_data, threads, timeout):
    if not isinstance(organism, str) or not organism.strip() or "\0" in organism:
        raise ValueError("organism must be a nonempty name")
    if receptor not in ("BCR", "TCR"):
        raise ValueError("receptor must be BCR or TCR")
    if not germline_v or not germline_j:
        raise ValueError("IgBLAST requires germline V and J database prefixes")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    if isinstance(timeout, bool) or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be finite and positive")
    if igdata and not (Path(igdata).expanduser() / "internal_data").is_dir():
        raise ValueError("IGDATA must contain the IgBLAST internal_data directory")
    _resolve_auxiliary(auxiliary_data, igdata)


def validate_annotation_options(*, executable="igblastn", organism="human",
                                receptor="BCR", germline_v, germline_j,
                                germline_d=None, auxiliary_data=None,
                                germline_c=None, igdata=None, threads=1, timeout=300):
    """Preflight annotation options and return the resolved executable path.

    Checks scalar parameters, executable availability, an explicitly supplied
    IGDATA directory, and readability of a supplied J auxiliary file.
    Explicit database paths must have discoverable nucleotide database files.
    Bare names may be resolved by BLASTDB or other IgBLAST configuration at run
    time. File presence does not assert that a database is valid or complete.
    """
    _validate_parameters(organism, receptor, germline_v, germline_j, igdata,
                         auxiliary_data, threads, timeout)
    for role, value in (("V", germline_v), ("D", germline_d), ("J", germline_j), ("C", germline_c)):
        if value and _explicit_prefix(value) and not _find_database(value)[1]:
            raise ValueError("IgBLAST %s database prefix has no nucleotide database files: %s" % (role, value))
    program = shutil.which(str(Path(executable).expanduser()))
    if not program:
        raise AnnotationError("Cannot find executable IgBLAST: " + str(executable))
    return str(Path(program).resolve())


def annotate_fasta(fasta_path, output_path, *, executable="igblastn", organism="human",
                   receptor="BCR", germline_v, germline_j, germline_d=None,
                   auxiliary_data=None, germline_c=None, igdata=None, threads=1,
                   timeout=300, records=None):
    """Run local IgBLAST using AIRR ``-outfmt 19`` and publish valid output.

    Germline arguments are BLAST database prefixes (not raw FASTA filenames).
    ``igdata`` sets IGDATA only in the child process. J auxiliary annotations
    must match the J database to obtain reliable frame and CDR3 information.
    If supplied, ``records`` contributes only cell_id and basic_* metadata,
    joined by sequence_id. Biological fields always come from IgBLAST.
    Existing output survives subprocess failure, timeout, or invalid output.
    """
    fasta_path, output_path = Path(fasta_path).resolve(), Path(output_path).resolve()
    if fasta_path == output_path:
        raise ValueError("Annotation output must differ from the input FASTA")
    identifiers = _fasta_ids(fasta_path)
    if identifiers:
        program = validate_annotation_options(
            executable=executable, organism=organism, receptor=receptor,
            germline_v=germline_v, germline_j=germline_j, germline_d=germline_d,
            germline_c=germline_c, auxiliary_data=auxiliary_data, igdata=igdata,
            threads=threads, timeout=timeout)
    else:
        # Empty biological output needs no installed executable, but malformed
        # options or metadata should not be silently accepted.
        _validate_parameters(organism, receptor, germline_v, germline_j, igdata,
                             auxiliary_data, threads, timeout)
    auxiliary_data = _resolve_auxiliary(auxiliary_data, igdata)
    assembly = {}
    for record in records or []:
        row = _normalize(record)
        identifier = row.get("sequence_id")
        if not identifier or identifier in assembly:
            raise ValueError("Assembly metadata requires unique, nonempty sequence_id values")
        assembly[identifier] = row
    if set(assembly) - set(identifiers):
        raise ValueError("Assembly metadata contains sequence IDs absent from input FASTA")
    metadata = {"tool": "igblastn", "format": "AIRR TSV", "organism": organism,
                "receptor": receptor, "germline_v": str(germline_v),
                "germline_d": str(germline_d) if germline_d else None,
                "germline_j": str(germline_j),
                "germline_c": str(germline_c) if germline_c else None,
                "auxiliary_data": str(auxiliary_data) if auxiliary_data else None,
                "igdata": str(igdata) if igdata else os.environ.get("IGDATA"),
                "working_directory": os.getcwd(), "warnings": []}
    if not auxiliary_data:
        metadata["warnings"].append("No J auxiliary data supplied; CDR3/frame annotations may be unavailable or incomplete.")
    if not germline_d:
        metadata["warnings"].append("No D database explicitly selected; supply a matching D database for VDJ loci.")
    if not identifiers:
        write_airr([], output_path)
        metadata.update(status="skipped_empty_fasta", record_count=0, elapsed_seconds=0)
        return metadata
    metadata["references"] = annotation_reference_fingerprints(
        germline_v=germline_v, germline_j=germline_j, germline_d=germline_d,
        germline_c=germline_c, auxiliary_data=auxiliary_data, igdata=igdata,
        organism=organism, executable=program)
    environment = os.environ.copy()
    if igdata:
        igdata_path = Path(igdata).expanduser().resolve()
        environment["IGDATA"] = str(igdata_path)
    started = time.monotonic()
    # Same filesystem as destination; no shell, shared temporary name, or append.
    with tempfile.TemporaryDirectory(prefix=".basic-igblast-", dir=output_path.parent) as directory:
        temporary_output = Path(directory) / "annotation.tsv"
        command = [program, "-query", str(fasta_path), "-out", str(temporary_output),
                   "-outfmt", "19", "-organism", str(organism), "-ig_seqtype",
                   "Ig" if receptor == "BCR" else "TCR", "-germline_db_V",
                   str(germline_v), "-germline_db_J", str(germline_j),
                   "-num_threads", str(threads)]
        for flag, value in (("-germline_db_D", germline_d),
                            ("-c_region_db", germline_c),
                            ("-auxiliary_data", auxiliary_data)):
            if value:
                command.extend((flag, str(value)))
        try:
            result = _run_igblast(command, environment, timeout)
        except subprocess.TimeoutExpired as error:
            raise AnnotationError("IgBLAST exceeded its %s second timeout" % timeout) from error
        except OSError as error:
            raise AnnotationError("Could not run IgBLAST: " + str(error)) from error
        if result.returncode:
            raise AnnotationError("IgBLAST exited with status %d: %s" %
                                  (result.returncode, result.stderr.strip()[-4000:]))
        if not temporary_output.is_file():
            raise AnnotationError("IgBLAST completed without producing AIRR output")
        annotated = read_airr(temporary_output)
        if {row["sequence_id"] for row in annotated} != set(identifiers):
            raise AnnotationError("IgBLAST AIRR sequence IDs do not match the input FASTA")
        for row in annotated:
            source = assembly.get(row["sequence_id"], {})
            for key, value in source.items():
                if key == "cell_id" or key.startswith("basic_"):
                    row[key] = value
            if source.get("locus"):
                row["basic_seed_locus"] = source["locus"]
            row["basic_locus_source"] = "igblast" if row.get("locus") else ""
        write_airr(annotated, output_path)
        metadata.update(status="complete", record_count=len(annotated), command=command,
                        stderr=result.stderr.strip(), elapsed_seconds=time.monotonic() - started)
    try:
        version_result = _run_igblast([program, "-version"], environment, min(timeout, 10))
        metadata["version"] = ((version_result.stdout.strip().splitlines() or [None])[0]
                               if version_result.returncode == 0 else None)
    except (OSError, subprocess.TimeoutExpired):
        metadata["version"] = None
    return metadata
