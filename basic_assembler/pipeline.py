"""Single-sample orchestration, provenance, and transactional output files."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from functools import lru_cache
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import tempfile
import time

from . import __version__
from .alignment import run_alignment
from .annotation import annotate_fasta, write_airr
from .assembly import AssemblyOptions, OverlapIndex, assemble_chain
from .reads import load_read_pool


class AssemblyIncomplete(RuntimeError):
    pass


def default_database_path() -> Path:
    checkout = Path(__file__).resolve().parent.parent / 'db'
    if (checkout / 'human').is_dir():
        return checkout
    from importlib.resources import files
    return Path(str(files('basic_assembler.references')))


@lru_cache(maxsize=512)
def _file_hash(path: str, size: int, modified: int) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(path) -> dict:
    path = Path(path).resolve()
    info = path.stat()
    return {'path': str(path), 'bytes': info.st_size,
            'sha256': _file_hash(str(path), info.st_size, info.st_mtime_ns)}


def atomic_json(path, value):
    path = Path(path)
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix=f'.{path.name}.', delete=False) as handle:
            tmp = Path(handle.name)
            json.dump(value, handle, indent=2, sort_keys=True, default=str, allow_nan=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


@contextmanager
def sample_lock(output: Path, name: str):
    path = output / f'.{name}.basic.lock'
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f'Another run owns {path}. If a previous process was killed, remove the stale lock after confirming it is no longer running.') from exc
    try:
        with os.fdopen(fd, 'w') as handle:
            handle.write(f'pid={os.getpid()}\n')
        yield
    finally:
        path.unlink(missing_ok=True)


def publish_bundle(stage: Path, output: Path, names: list[str]):
    """Replace a complete bundle, rolling back on a normal publication error.

    Each rename is atomic; run.json is published last as the completion marker.
    Readers should verify the output hashes in that manifest. A power loss can
    interrupt a multi-file transaction, so this is not a database transaction.
    """
    # Outside stage: its TemporaryDirectory must not delete originals if the
    # filesystem also refuses a rollback operation.
    backup = Path(tempfile.mkdtemp(prefix='.basic-recovery-', dir=output))
    installed, backed_up = [], []
    try:
        for name in names:
            target = output / name
            if target.exists():
                os.replace(target, backup / name)
                backed_up.append(name)
            os.replace(stage / name, target)
            installed.append(name)
    except BaseException as original_error:
        recovery_errors = []
        for name in reversed(installed):
            target = output / name
            try:
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink(missing_ok=True)
            except OSError as error:
                recovery_errors.append(str(error))
        for name in reversed(backed_up):
            try:
                os.replace(backup / name, output / name)
            except OSError as error:
                recovery_errors.append(str(error))
        if recovery_errors:
            raise RuntimeError(f'Publication and rollback failed. Previous files remain recoverable in {backup}: '
                               + '; '.join(recovery_errors)) from original_error
        shutil.rmtree(backup, ignore_errors=True)
        raise
    else:
        shutil.rmtree(backup, ignore_errors=True)


def _aggregate(candidates: list[dict]) -> list[dict]:
    """Count identical seed sequences across alleles of the same locus."""
    groups = {}
    for candidate in candidates:
        key = candidate['sequence'], candidate.get('locus', '')
        if key not in groups:
            groups[key] = dict(candidate, count=0, references=[])
        group = groups[key]
        group['count'] += candidate['count']
        group['mapq'] = max(group.get('mapq', 0), candidate.get('mapq', 0))
        group['references'].append(candidate['reference'])
    return sorted(groups.values(), key=lambda c: (-c['count'], -c['mapq'], c['sequence'], c.get('locus', '')))


def _seed_pairs(variable: list[dict], constant: list[dict], limit: int) -> list[tuple]:
    def selected(candidates):
        counts, result = {}, []
        for candidate in _aggregate(candidates):
            locus = candidate.get('locus', '')
            if counts.get(locus, 0) < limit:
                result.append(candidate)
                counts[locus] = counts.get(locus, 0) + 1
        return result

    starts, ends = selected(variable), selected(constant)
    if not starts:
        starts = [None]
    if not ends:
        ends = [None]
    pairs = [(v, c) for v in starts for c in ends
             if v is None or c is None or not v.get('locus') or not c.get('locus')
             or v['locus'] == c['locus']]
    # Never join incompatible loci; retain unmatched one-sided evidence too.
    paired_v = {id(v) for v, c in pairs}
    paired_c = {id(c) for v, c in pairs}
    pairs.extend((v, None) for v in starts if v is not None and id(v) not in paired_v)
    pairs.extend((None, c) for c in ends if c is not None and id(c) not in paired_c)
    pairs.sort(key=lambda pair: (-sum(x['count'] for x in pair if x),
                                 tuple(x['sequence'] if x else '' for x in pair)))
    # Selection already bounds each locus to ``limit`` candidates per side.
    # Keep every selected locus, including four-locus custom TCR references.
    # At most len(starts)*len(ends) + len(starts) + len(ends) pairs survive.
    return pairs


def _assemble_results(args, pool, mappings):
    options = AssemblyOptions(**{key: args[key] for key in
                              ('min_overlap', 'min_support', 'max_steps', 'max_contig_length', 'scoring', 'beam_width')})
    index = OverlapIndex(pool.counts)
    chain_names = ('heavy', 'light') if args['type'] == 'BCR' else _tcr_chain_names(args, mappings)
    chains, records = [], []
    for name, vregion, cregion in zip(chain_names, ('hv', 'lv'), ('hc', 'lc')):
        attempts = []
        for variable, constant in _seed_pairs(mappings[vregion].candidates, mappings[cregion].candidates, args['max_anchors']):
            result = assemble_chain(variable['sequence'] if variable else '',
                                    constant['sequence'] if constant else '',
                                    pool.counts, options, index=index,
                                    quality_weights=pool.quality_weights if args['quality_weighted'] else None,
                                    mate_links=pool.mate_links if args['mate_aware'] else None)
            actual_name = ((variable or constant or {}).get('locus') or name) if args['type'] == 'TCR' else name
            result.update(chain=actual_name, variable_anchor=variable, constant_anchor=constant)
            attempts.append(result)
        attempts.sort(key=lambda a: (a['status'] != 'complete',
                                     -sum(x['count'] for x in (a['variable_anchor'], a['constant_anchor']) if x),
                                     tuple(c['sequence'] for c in a['contigs'])))
        selected, seen = [], set()
        for attempt in attempts:
            if not args['allow_partial'] and attempts[0]['status'] == 'complete' and attempt['status'] != 'complete':
                continue
            signature = tuple((c['region'], c['sequence']) for c in attempt['contigs'])
            if signature in seen:
                continue
            seen.add(signature)
            selected.append(attempt)
            if len(selected) >= args['max_candidates']:
                break
        best_status = selected[0]['status'] if selected else 'unsupported'
        best_name = selected[0]['chain'] if selected else name
        chains.append({'chain': best_name, 'region_group': vregion, 'status': best_status, 'candidates': selected,
                       'attempted_seed_pairs': len(attempts)})
        for number, result in enumerate(selected, 1):
            v, c = result['variable_anchor'], result['constant_anchor']
            locus = (v or c or {}).get('locus', '')
            region_counts = {}
            for contig in result['contigs']:
                region_counts[contig['region']] = region_counts.get(contig['region'], 0) + 1
                suffix = '' if contig['region'] == 'full' else f"[{contig['region']}_region_contig]"
                identifier = f"cell_id={args['name']};{result['chain']}_chain{suffix};BASIC"
                if number > 1:
                    identifier += f';candidate={number}'
                if region_counts[contig['region']] > 1:
                    identifier += f";path={region_counts[contig['region']]}"
                if any(r['sequence_id'] == identifier for r in records):
                    identifier += f';group={vregion}'
                records.append({'sequence_id': identifier, 'sequence': contig['sequence'],
                                'cell_id': args['name'], 'locus': locus,
                                'status': result['status'], 'ambiguous': contig.get('ambiguous', False),
                                'basic_chain': result['chain'], 'basic_region_group': vregion, 'basic_region': contig['region'],
                                'basic_target_species': args['genome'],
                                'basic_seed_species': args['seed_species'] or args['genome'],
                                'basic_reference_transfer': bool(args['seed_species'] and args['seed_species'] != args['genome']),
                                'basic_candidate': number, 'basic_locus_source': 'seed_alignment' if locus else '',
                                'basic_anchor_support': min(x['count'] for x in (v, c) if x),
                                'basic_path_min_support': contig.get('read_support_min'),
                                'basic_termination': contig.get('termination', '')})
    return chains, records


def _tcr_chain_names(args, mappings):
    """The historical mouse and human databases use opposite h/l order."""
    names = []
    for number, region in enumerate(('hv', 'lv'), 1):
        counterpart = 'hc' if region == 'hv' else 'lc'
        pairs = _seed_pairs(mappings[region].candidates, mappings[counterpart].candidates, 1)
        observed = [(v or c or {}).get('locus') for v, c in pairs]
        if observed and observed[0]:
            names.append(observed[0])
            continue
        loci = {c.get('locus') for c in mappings[region].candidates if c.get('locus')}
        reference = Path(args['db_path']) / region
        if reference.is_file():
            with reference.open(encoding='utf-8') as handle:
                reference_loci = {match.group(1) for line in handle if line.startswith('>')
                                  for match in re.finditer(r'(TR[ABGD])[VCJD]', line)}
            if reference_loci:
                loci = reference_loci
        names.append(next(iter(loci)) if len(loci) == 1 else f'TCR_{number}')
    if names[0] == names[1]:
        return ('TCR_1', 'TCR_2')
    return tuple(names)


def _summary(path, chains, records):
    fields = ['chain', 'status', 'contigs', 'ambiguous_contigs', 'attempted_seed_pairs']
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for chain in chains:
            matching = [r for r in records if r['basic_region_group'] == chain['region_group']]
            writer.writerow({'chain': chain['chain'], 'status': chain['status'],
                             'contigs': len(matching),
                             'ambiguous_contigs': sum(bool(r['ambiguous']) for r in matching),
                             'attempted_seed_pairs': chain['attempted_seed_pairs']})


def run_sample(args: dict) -> dict:
    output = Path(args['output_location'])
    output.mkdir(parents=True, exist_ok=True)
    name = args['name']
    filenames = [f'{name}.fasta', f'{name}.airr.tsv', f'{name}.summary.tsv', f'{name}.logs', f'{name}.run.json']
    with sample_lock(output, name):
        existing = [p for p in filenames if (output / p).exists()]
        if existing and not args['overwrite']:
            raise FileExistsError(f"Output already exists: {', '.join(existing)}. Choose a new -n or use --overwrite.")
        if args['tmpdir']:
            Path(args['tmpdir']).mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix=f'basic-{name}-', dir=args['tmpdir']))
        started = time.monotonic()
        manifest = {'basic_version': __version__, 'sample_id': name, 'status': 'running',
                    'started_utc': datetime.now(timezone.utc).isoformat(),
                    'python_version': platform.python_version(), 'platform': platform.platform(),
                    'parameters': args, 'bowtie_version': args['bowtie_version'],
                    'experimental': args['scoring'] != 'legacy' or args['beam_width'] > 1
                                    or args['max_anchors'] > 1 or args['max_candidates'] > 1
                                    or args['local_seeding'] or bool(args['seed_species'] and args['seed_species'] != args['genome']),
                    'target_species': args['genome'], 'seed_species': args['seed_species'] or args['genome']}
        try:
            manifest['inputs'] = [fingerprint(p) for p in args['single'] + args['left'] + args['right']]
            manifest['references'] = [fingerprint(p) for p in sorted(Path(args['db_path']).iterdir())
                                      if p.is_file() and (p.suffix in ('.bt2', '.bt2l') or p.name in ('hv', 'hc', 'lv', 'lc', 'reference_manifest.json'))]
            prepared = work / 'reads.fastq'
            pool = load_read_pool(args['single'], args['left'], args['right'],
                                  subsample=args['subsample'], min_read_length=args['min_overlap'] + 1,
                                  trim_quality=args['trim_quality'], quality_weighted=args['quality_weighted'],
                                  prepared_path=prepared, retain_mates=args['mate_aware'])
            manifest['read_statistics'] = pool.stats
            manifest['prepared_input'] = fingerprint(prepared)
            if args['VERBOSE']:
                print(f"{name}: retained {pool.stats['retained_reads']} reads; mapping four reference regions", flush=True)
            mappings = {region: run_alignment(args['bowtie'], args['indexes'][region], [str(prepared)],
                                              args['num_threads'], work / f'{region}.stderr.log', strand=region,
                                              local=args['local_seeding'], timeout=args['alignment_timeout'],
                                              original_sequences=pool.ambiguous_sequences)
                        for region in ('hv', 'hc', 'lv', 'lc')}
            manifest['mapping'] = {region: {'mapped_reads': result.mapped_reads, 'command': result.command,
                                           'candidate_seeds': len(result.candidates)}
                                   for region, result in mappings.items()}
            chains, records = _assemble_results(args, pool, mappings)
            manifest['chains'] = chains
            status = 'complete' if all(c['status'] == 'complete' for c in chains) else ('partial' if records else 'unsupported')
            manifest['status'] = status
            if status != 'complete' and not args['allow_partial']:
                detail = ', '.join(f"{c['chain']}={c['status']}" for c in chains)
                raise AssemblyIncomplete(f'Assembly incomplete ({detail}); use -a to retain partial results')
            with tempfile.TemporaryDirectory(prefix=f'.{name}.publish-', dir=output) as stage_dir:
                stage = Path(stage_dir)
                fasta = stage / f'{name}.fasta'
                with fasta.open('w', encoding='utf-8') as handle:
                    for record in records:
                        handle.write(f">{record['sequence_id']}\n{record['sequence']}\n")
                airr = stage / f'{name}.airr.tsv'
                if args['igblast'] and records:
                    manifest['annotation'] = annotate_fasta(
                        fasta, airr, executable=args['igblast'], organism=args['genome'], receptor=args['type'],
                        germline_v=args['germline_v'], germline_d=args['germline_d'], germline_j=args['germline_j'],
                        germline_c=args['germline_c'], auxiliary_data=args['auxiliary_data'], igdata=args['igdata'],
                        threads=args['num_threads'], timeout=args['annotation_timeout'], records=records)
                else:
                    write_airr(records, airr)
                    manifest['annotation'] = {'status': 'not_run', 'reason': 'no_contigs' if args['igblast'] else 'not_requested'}
                _summary(stage / f'{name}.summary.tsv', chains, records)
                logs = stage / f'{name}.logs'
                logs.mkdir()
                for log in work.glob('*.log'):
                    shutil.copy2(log, logs / log.name)
                manifest['elapsed_seconds'] = time.monotonic() - started
                manifest['outputs'] = [{**fingerprint(stage / filename), 'path': str(output / filename)}
                                       for filename in filenames[:-1] if (stage / filename).is_file()]
                if args['keep_temp']:
                    manifest['temporary_directory'] = str(work)
                atomic_json(stage / f'{name}.run.json', manifest)
                publish_bundle(stage, output, filenames)
            (output / f'{name}.failed.json').unlink(missing_ok=True)
            return {'sample_id': name, 'status': status, 'fasta': str(output / f'{name}.fasta'),
                    'manifest': str(output / f'{name}.run.json'), 'elapsed_seconds': manifest['elapsed_seconds']}
        except BaseException as exc:
            manifest.update(status='failed', error=str(exc), error_type=type(exc).__name__,
                            elapsed_seconds=time.monotonic() - started)
            # Preserve diagnostics, while leaving a previous successful bundle intact.
            manifest['stderr'] = {p.name: p.read_text(errors='replace')[-16000:] for p in work.glob('*.log')}
            if args['keep_temp']:
                manifest['temporary_directory'] = str(work)
            atomic_json(output / f'{name}.failed.json', manifest)
            raise
        finally:
            if not args['keep_temp']:
                shutil.rmtree(work)


def write_batch_summary(results: list[dict], output_location):
    output = Path(output_location).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    # Unique summary paths avoid collisions between concurrent batch invocations.
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    path = output / f'batch-{stamp}.json'
    atomic_json(path, {'basic_version': __version__, 'samples': results,
                       'failed': sum(r['status'] == 'failed' for r in results)})
    print(f'Batch summary: {path}')
