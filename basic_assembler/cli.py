"""Command-line interface; the original BASIC.py flags remain available."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path

from . import __version__


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Assemble BCR/TCR sequences from per-cell FASTQ reads.")
    p.add_argument('-i', dest='type', choices=['BCR', 'TCR'], default='BCR')
    p.add_argument('-p', dest='num_threads', type=int, default=2, help='Bowtie2 threads (default: 2)')
    p.add_argument('-n', dest='name', default='result', help='Sample/output basename')
    p.add_argument('-SE', dest='FASTQ', default='', help='Single-end FASTQ file(s), comma-separated')
    p.add_argument('-PE_1', dest='LEFT', default='', help='Mate 1 FASTQ file(s), comma-separated')
    p.add_argument('-PE_2', dest='RIGHT', default='', help='Mate 2 FASTQ file(s), comma-separated')
    p.add_argument('-g', dest='genome', default='human', help='Species directory under the reference root')
    p.add_argument('--seed-species', help='Use another species reference directory for exploratory seeding')
    p.add_argument('--local-seeding', action='store_true', help='Experimental sensitive local homology search for seed reads')
    p.add_argument('-b', dest='bowtie', default='bowtie2', help='Bowtie2 executable or directory; searches PATH')
    p.add_argument('-t', dest='tmpdir', default=None, help='Parent for unique temporary run directories')
    p.add_argument('-o', dest='output_location', default='.', help='Output directory')
    p.add_argument('-a', dest='allow_partial', action='store_true', help='Write partial or empty biological results')
    p.add_argument('-v', dest='VERBOSE', action='store_true')
    p.add_argument('-d', dest='database_path', default='', help='Reference root containing species/receptor directories')
    p.add_argument('--subsample', type=int, help='Use the first N reads, or N synchronized pairs, across input lanes')
    p.add_argument('--overwrite', action='store_true', help='Replace an existing successful output bundle')
    p.add_argument('--keep-temp', action='store_true', help='Keep prepared reads and Bowtie2 logs')
    p.add_argument('--check', action='store_true', help='Validate inputs and dependencies without assembling')
    p.add_argument('--batch', type=Path, help='TSV with sample_id and fastq, or left/right; paths relative to TSV')
    p.add_argument('--continue-on-error', action='store_true', help='In a batch, finish other samples after failures')
    p.add_argument('--min-overlap', type=int, default=15)
    p.add_argument('--min-support', type=int, default=2)
    p.add_argument('--max-steps', type=int, default=10000)
    p.add_argument('--max-contig-length', type=int, default=30000)
    p.add_argument('--scoring', choices=['legacy', 'evidence'], default='legacy', help='Evidence scoring is experimental')
    p.add_argument('--beam-width', type=int, default=1, help='Bounded alternative path search; >1 is experimental')
    p.add_argument('--max-anchors', type=int, default=1, help='Candidate seeds per region/locus; >1 is experimental')
    p.add_argument('--max-candidates', type=int, default=1, help='Maximum distinct assemblies per chain')
    p.add_argument('--trim-quality', type=int, default=0, help='Trim terminal Phred+33 bases below this quality')
    p.add_argument('--quality-weighted', action='store_true', help='Use read reliability in evidence scoring')
    p.add_argument('--mate-aware', action='store_true', help='Use mate support in evidence scoring')
    p.add_argument('--igblast', metavar='EXECUTABLE', help='Optionally annotate assembled FASTA with igblastn')
    p.add_argument('--germline-v', help='IgBLAST V database prefix')
    p.add_argument('--germline-d', help='IgBLAST D database prefix')
    p.add_argument('--germline-j', help='IgBLAST J database prefix')
    p.add_argument('--germline-c', help='Optional IgBLAST C database prefix')
    p.add_argument('--auxiliary-data', help='Species-matched IgBLAST J auxiliary data')
    p.add_argument('--igdata', help='IgBLAST data directory containing internal_data')
    p.add_argument('--annotation-timeout', type=int, default=300)
    p.add_argument('--alignment-timeout', type=float, help='Optional per-region mapping timeout in seconds')
    p.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return p


def parse_args(argv=None) -> dict:
    return vars(parser().parse_args(argv))


def _input_paths(value: str, base: Path | None = None) -> list[str]:
    if not value:
        return []
    paths = value.split(',')
    if any(not s.strip() for s in paths):
        raise ValueError('Empty filename in a comma-separated input list')
    result = []
    for name in paths:
        path = Path(name).expanduser()
        if base is not None and not path.is_absolute():
            path = base / path
        if not path.is_file():
            raise ValueError(f'Input file not found: {path}')
        result.append(str(path.resolve()))
    return result


def _output_overlaps_input(source: Path, destination: Path) -> bool:
    if source == destination or destination in source.parents:
        return True
    # resolve() follows symlinks but need not canonicalize filename case.
    # File identity also covers differently cased parent directories that a
    # published logs directory would replace. Reject hardlinks conservatively.
    try:
        target_stat = destination.stat()
    except (FileNotFoundError, NotADirectoryError):
        return False
    for path in (source, *source.parents):
        try:
            if os.path.samestat(target_stat, path.stat()):
                return True
        except (FileNotFoundError, NotADirectoryError):
            continue
    return False


def config_args(values: dict) -> dict:
    """Validate an individual sample without creating output files."""
    from .alignment import resolve_bowtie, validate_indexes
    from .pipeline import default_database_path

    args = parse_args([])
    args.update(values)
    if args['type'] not in ('BCR', 'TCR'):
        raise ValueError('Receptor must be BCR or TCR')
    for key in ('name', 'genome', 'seed_species'):
        if key == 'seed_species' and not args[key]:
            continue
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args[key]):
            raise ValueError(f'{key} must be a basename using letters, numbers, dot, underscore or hyphen')
    for key in ('num_threads', 'min_overlap', 'min_support', 'max_steps', 'max_contig_length', 'beam_width', 'max_anchors', 'max_candidates', 'annotation_timeout'):
        if args[key] < 1:
            raise ValueError(f'{key} must be positive')
    if args['subsample'] is not None and args['subsample'] < 1:
        raise ValueError('subsample must be positive')
    if args['alignment_timeout'] is not None and (not math.isfinite(args['alignment_timeout']) or args['alignment_timeout'] <= 0):
        raise ValueError('alignment_timeout must be finite and positive')
    if args['max_contig_length'] < args['min_overlap']:
        raise ValueError('max_contig_length must be at least min_overlap')
    if not 0 <= args['trim_quality'] <= 93:
        raise ValueError('trim_quality must be between 0 and 93')
    if args['scoring'] != 'evidence' and (args['quality_weighted'] or args['mate_aware']):
        raise ValueError('--quality-weighted and --mate-aware require --scoring evidence')
    if args['FASTQ'] and (args['LEFT'] or args['RIGHT']):
        raise ValueError('Use either -SE or both -PE_1 and -PE_2, never a mixture')
    if not args['FASTQ'] and not (args['LEFT'] and args['RIGHT']):
        raise ValueError('Supply -SE or both -PE_1 and -PE_2')
    args['single'] = _input_paths(args['FASTQ'])
    args['left'] = _input_paths(args['LEFT'])
    args['right'] = _input_paths(args['RIGHT'])
    if len(args['left']) != len(args['right']):
        raise ValueError('Mate file lists must contain the same number of lanes')
    args['read_type'] = 'single' if args['single'] else 'paired'
    if args['mate_aware'] and args['read_type'] != 'paired':
        raise ValueError('--mate-aware requires paired input')
    if args['igblast'] and not (args['germline_v'] and args['germline_j']):
        raise ValueError('--igblast requires --germline-v and --germline-j')
    if not args['igblast'] and any(args[k] for k in ('germline_v', 'germline_d', 'germline_j', 'germline_c', 'auxiliary_data', 'igdata')):
        raise ValueError('Annotation options require --igblast')
    if args['igblast']:
        from .annotation import validate_annotation_options
        args['igblast'] = validate_annotation_options(
            executable=args['igblast'], organism=args['genome'], receptor=args['type'],
            germline_v=args['germline_v'], germline_d=args['germline_d'], germline_j=args['germline_j'],
            germline_c=args['germline_c'], auxiliary_data=args['auxiliary_data'], igdata=args['igdata'],
            threads=args['num_threads'], timeout=args['annotation_timeout'])
    args['database_path'] = str(Path(args['database_path'] or default_database_path()).expanduser().resolve())
    args['db_path'] = str(Path(args['database_path']) / (args['seed_species'] or args['genome']) / args['type'])
    args['indexes'] = validate_indexes(args['db_path'])
    args['bowtie'], args['bowtie_version'] = resolve_bowtie(args['bowtie'])
    args['output_location'] = str(Path(args['output_location']).expanduser().resolve())
    output = Path(args['output_location'])
    destinations = [output / f"{args['name']}{suffix}" for suffix in
                    ('.fasta', '.airr.tsv', '.summary.tsv', '.logs', '.run.json', '.failed.json')]
    for filename in args['single'] + args['left'] + args['right']:
        source = Path(filename).resolve()
        if any(_output_overlaps_input(source, destination.resolve()) for destination in destinations):
            raise ValueError(f'Output destination overlaps an input file: {source}')
    for key in ('output_location', 'tmpdir'):
        if args[key]:
            path = Path(args[key]).expanduser().resolve()
            if path.exists() and not path.is_dir():
                raise ValueError(f'{key} is not a directory: {path}')
            args[key] = str(path)
    return args


def read_batch(path: Path, defaults: dict) -> list[dict]:
    """Read and validate all identifiers before the first sample is run."""
    path = path.expanduser().resolve()
    samples, seen = [], set()
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or 'sample_id' not in reader.fieldnames:
            raise ValueError('Batch TSV requires a sample_id column')
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('Batch TSV has duplicate column names')
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'Incorrect number of TSV fields on batch line {number}')
            name = (row.get('sample_id') or '').strip()
            if not name or name.casefold() in seen:
                raise ValueError(f'Empty or duplicate sample_id on batch line {number}: {name!r}')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
                raise ValueError(f'Invalid sample_id on batch line {number}: {name!r}')
            seen.add(name.casefold())
            args = dict(defaults, name=name, batch=None)
            for column, key in (('fastq', 'FASTQ'), ('left', 'LEFT'), ('right', 'RIGHT')):
                value = (row.get(column) or '').strip()
                paths = []
                for item in value.split(',') if value else []:
                    if not item:
                        raise ValueError(f'Empty filename on batch line {number}')
                    candidate = Path(item).expanduser()
                    paths.append(str(candidate if candidate.is_absolute() else path.parent / candidate))
                args[key] = ','.join(paths)
            args['genome'] = row.get('genome') or defaults['genome']
            args['type'] = row.get('receptor') or defaults['type']
            samples.append(args)
    if not samples:
        raise ValueError('Batch TSV contains no samples')
    return samples


def _validate_batch_destinations(samples: list[dict], manifest: Path) -> None:
    """Reject cross-sample input collisions before publishing any batch output.

    Resolve paths without requiring them to exist: individual missing inputs still
    follow --continue-on-error. Keep lexical paths too, since replacing a logs
    directory can remove an input symlink even when its target is elsewhere.
    """
    sources = {manifest.expanduser().absolute(), manifest.expanduser().resolve()}
    for sample in samples:
        for key in ('FASTQ', 'LEFT', 'RIGHT'):
            for filename in (sample.get(key) or '').split(','):
                if filename.strip():
                    path = Path(filename).expanduser()
                    sources.update((path.absolute(), path.resolve()))
    for sample in samples:
        output = Path(sample['output_location']).expanduser()
        for suffix in ('.fasta', '.airr.tsv', '.summary.tsv', '.logs', '.run.json', '.failed.json'):
            destination = output / f"{sample['name']}{suffix}"
            for target in (destination.absolute(), destination.resolve()):
                if any(_output_overlaps_input(source, target) for source in sources):
                    raise ValueError(f'Batch output destination overlaps an input file or batch manifest: {destination}')


def main(argv=None) -> int:
    from .pipeline import run_sample, write_batch_summary

    args = parse_args(argv)
    try:
        if args['batch'] and (args['FASTQ'] or args['LEFT'] or args['RIGHT']):
            raise ValueError('--batch cannot be combined with -SE/-PE_1/-PE_2')
        samples = read_batch(args['batch'], args) if args['batch'] else [args]
        if args['batch']:
            _validate_batch_destinations(samples, args['batch'])
        results = []
        for sample in samples:
            try:
                config = config_args(sample)
                if args['check']:
                    result = {'sample_id': config['name'], 'status': 'checked', 'bowtie': config['bowtie'], 'bowtie_version': config['bowtie_version'], 'db_path': config['db_path']}
                    print(json.dumps(result, sort_keys=True))
                else:
                    result = run_sample(config)
                    print(f"{config['name']}: {result['status']} — {result['fasta']}")
                results.append(result)
            except (OSError, ValueError, RuntimeError) as exc:
                print(f"BASIC: {sample['name']}: {exc}", file=sys.stderr)
                results.append({'sample_id': sample['name'], 'status': 'failed', 'error': str(exc)})
                if not args['continue_on_error']:
                    break
        if args['batch'] and not args['check']:
            write_batch_summary(results, args['output_location'])
        return 1 if any(r['status'] == 'failed' for r in results) else 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'BASIC: {exc}', file=sys.stderr)
        return 1
