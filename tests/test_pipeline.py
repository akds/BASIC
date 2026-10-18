import csv
import json
import random
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from basic_assembler import pipeline
from basic_assembler.annotation import read_airr
from basic_assembler.assembly import revcom
from basic_assembler.cli import config_args, main, parse_args, read_batch


@pytest.fixture
def sample(tmp_path):
    rng = random.Random(872)
    targets = []
    while len(targets) < 2:
        target = ''.join(rng.choices('ACGT', k=96))
        if all(max(target[i:i + 24].count(base) for base in 'ACGT') < 12 for i in range(0, 73, 8)):
            targets.append(target)
    reads = tmp_path / 'input with spaces +ü.fastq'
    with reads.open('w') as handle:
        i = 0
        for target in targets:
            for offset in range(0, 73, 8):
                for _ in range(3):
                    seq = target[offset:offset + 24]
                    handle.write(f'@read{i}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
                    i += 1
    seeds = {'hv': targets[0][:24], 'hc': revcom(targets[0][-24:]),
             'lv': targets[1][:24], 'lc': revcom(targets[1][-24:])}
    program = tmp_path / 'bin with spaces' / 'bowtie2'
    program.parent.mkdir()
    program.write_text(f'''#!{sys.executable}
import sys
from pathlib import Path
if '--version' in sys.argv:
    print('bowtie2 fixture 1')
    raise SystemExit(0)
seeds = {seeds!r}
region = Path(sys.argv[sys.argv.index('-x')+1]).name
seq = seeds[region]
locus = 'IGH' if region.startswith('h') else 'IGK'
for i in range(3):
    print(f'read{{i}}\t0\t{{locus}}V1\t1\t42\t24M\t*\t0\t0\t{{seq}}\t'+'I'*24)
''')
    program.chmod(0o755)
    indexes = tmp_path / 'db' / 'human' / 'BCR'
    indexes.mkdir(parents=True)
    for region in seeds:
        for suffix in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
            (indexes / f'{region}.{suffix}.bt2').write_text('fixture')
    args = ['-SE', str(reads), '-b', str(program), '-d', str(indexes.parent.parent),
            '-o', str(tmp_path / 'output'), '-t', str(tmp_path / 'temp')]
    return SimpleNamespace(args=args, reads=reads, program=program, targets=targets,
                           output=tmp_path / 'output', tmp=tmp_path)


def test_complete_output_bundle_is_consistent_and_reproducible(sample):
    assert main(sample.args) == 0
    manifest = json.loads((sample.output / 'result.run.json').read_text())
    assert manifest['status'] == 'complete'
    assert manifest['read_statistics']['consumed_reads'] == 60
    records = read_airr(sample.output / 'result.airr.tsv')
    assert len(records) == 2
    assert {r['sequence'] for r in records} == set(sample.targets)
    assert all(r['productive'] == '' and r['junction'] == '' for r in records)
    for record in manifest['outputs']:
        assert pipeline.fingerprint(record['path'])['sha256'] == record['sha256']
    original = (sample.output / 'result.fasta').read_bytes()
    assert main(sample.args) == 1
    assert (sample.output / 'result.fasta').read_bytes() == original
    assert main(sample.args + ['--overwrite']) == 0
    assert (sample.output / 'result.fasta').read_bytes() == original
    assert not list((sample.tmp / 'temp').iterdir())


def test_aligner_failure_is_not_partial_success(sample):
    sample.program.write_text(f'#!{sys.executable}\nimport sys\nif "--version" in sys.argv:\n print("fixture 1")\nelse:\n print("alignment broke",file=sys.stderr)\n sys.exit(7)\n')
    assert main(sample.args + ['-a']) == 1
    assert not (sample.output / 'result.fasta').exists()
    failure = json.loads((sample.output / 'result.failed.json').read_text())
    assert failure['status'] == 'failed'
    assert 'exit 7' in failure['error']
    assert 'alignment broke' in failure['stderr']['hv.stderr.log']
    assert not (sample.output / '.result.basic.lock').exists()


def test_failed_overwrite_preserves_previous_success(sample):
    assert main(sample.args) == 0
    previous = (sample.output / 'result.fasta').read_bytes()
    sample.reads.write_text('@bad\nACGT\n+\nI\n')
    assert main(sample.args + ['--overwrite']) == 1
    assert (sample.output / 'result.fasta').read_bytes() == previous
    assert (sample.output / 'result.failed.json').exists()


def test_disconnected_anchors_require_partial_option(sample, monkeypatch):
    def disconnected(*args, **kwargs):
        return {'status': 'partial', 'contigs': [{'sequence': 'ACGT' * 8, 'region': 'variable', 'ambiguous': False}], 'diagnostics': {}}
    monkeypatch.setattr(pipeline, 'assemble_chain', disconnected)
    assert main(sample.args) == 1
    assert not (sample.output / 'result.fasta').exists()
    assert main(sample.args + ['-a']) == 0
    assert json.loads((sample.output / 'result.run.json').read_text())['status'] == 'partial'


def test_multiple_beam_paths_have_unique_identifiers(sample, monkeypatch):
    def two_paths(*args, **kwargs):
        return {'status': 'complete', 'contigs': [
            {'sequence': 'ACGT' * 8, 'region': 'full', 'ambiguous': True},
            {'sequence': 'TGCA' * 8, 'region': 'full', 'ambiguous': True}], 'diagnostics': {}}
    monkeypatch.setattr(pipeline, 'assemble_chain', two_paths)
    assert main(sample.args + ['--beam-width', '2']) == 0
    records = read_airr(sample.output / 'result.airr.tsv')
    assert len({r['sequence_id'] for r in records}) == 4


def test_batch_continue_on_missing_file_runs_valid_sample(sample):
    batch = sample.tmp / 'cells.tsv'
    with batch.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample_id', 'fastq'])
        writer.writerow(['missing', 'nonexistent.fastq'])
        writer.writerow(['good', sample.reads.name])
    args = sample.args[2:] + ['--batch', str(batch), '--continue-on-error']
    assert main(args) == 1
    assert (sample.output / 'good.fasta').exists()
    summary = json.loads(next(sample.output.glob('batch-*.json')).read_text())
    assert [r['status'] for r in summary['samples']] == ['failed', 'complete']


@pytest.mark.parametrize('relative_input', ['first.fasta', 'first.failed.json', 'first.logs/reads.fastq'])
def test_batch_rejects_outputs_overlapping_later_input_before_running(sample, relative_input):
    later_input = sample.output / relative_input
    later_input.parent.mkdir(parents=True, exist_ok=True)
    original = sample.reads.read_bytes()
    later_input.write_bytes(original)
    batch = sample.tmp / 'cells.tsv'
    batch.write_text(f'sample_id\tfastq\nfirst\t{sample.reads}\nsecond\t{later_input}\n')
    args = sample.args[2:] + ['--batch', str(batch), '--overwrite', '--continue-on-error']
    assert main(args) == 1
    assert later_input.read_bytes() == original
    assert not (sample.output / 'first.run.json').exists()
    assert not (sample.output / 'second.run.json').exists()


def test_batch_preserves_input_symlink_inside_replaced_logs_directory(sample):
    logs = sample.output / 'first.logs'
    logs.mkdir(parents=True)
    later_input = logs / 'reads.fastq'
    later_input.symlink_to(sample.reads)
    batch = sample.tmp / 'cells.tsv'
    batch.write_text(f'sample_id\tfastq\nfirst\t{sample.reads}\nsecond\t{later_input}\n')
    assert main(sample.args[2:] + ['--batch', str(batch), '--overwrite']) == 1
    assert later_input.is_symlink()
    assert later_input.read_bytes() == sample.reads.read_bytes()
    assert not (sample.output / 'first.run.json').exists()


def test_batch_manifest_cannot_be_replaced_by_sample_output(sample):
    sample.output.mkdir()
    batch = sample.output / 'first.run.json'
    original = f'sample_id\tfastq\nfirst\t{sample.reads}\n'
    batch.write_text(original)
    assert main(sample.args[2:] + ['--batch', str(batch), '--overwrite']) == 1
    assert batch.read_text() == original
    assert not (sample.output / 'first.fasta').exists()


def test_batch_conservatively_rejects_input_hardlink_to_an_output_path(tmp_path):
    from basic_assembler.cli import _validate_batch_destinations
    source = tmp_path / 'reads.fastq'
    source.write_text('@read\nACGT\n+\nIIII\n')
    (tmp_path / 'first.fasta').hardlink_to(source)
    # Identity checks also protect aliases on case-insensitive filesystems.
    with pytest.raises(ValueError, match='overlaps'):
        _validate_batch_destinations([dict(name='first', output_location=str(tmp_path), FASTQ=str(source))],
                                     tmp_path / 'cells.tsv')
    assert source.read_text().startswith('@read')


@pytest.mark.parametrize('relative_input', ['first.FASTA', 'first.LOGS/reads.fastq'])
def test_batch_rejects_case_alias_of_later_input(sample, relative_input):
    later_input = sample.output / relative_input
    later_input.parent.mkdir(parents=True, exist_ok=True)
    original = sample.reads.read_bytes()
    later_input.write_bytes(original)
    if not (sample.output / relative_input.lower()).exists():
        pytest.skip('Filesystem is case-sensitive')
    batch = sample.tmp / 'cells.tsv'
    batch.write_text(f'sample_id\tfastq\nfirst\t{sample.reads}\nsecond\t{later_input}\n')
    assert main(sample.args[2:] + ['--batch', str(batch), '--overwrite']) == 1
    assert later_input.read_bytes() == original
    assert not (sample.output / 'first.run.json').exists()


def test_single_sample_rejects_case_alias_of_input(sample):
    original = sample.reads.read_bytes()
    sample.output.mkdir()
    source = sample.output / 'result.FASTA'
    source.write_bytes(original)
    if not (sample.output / 'result.fasta').exists():
        pytest.skip('Filesystem is case-sensitive')
    args = sample.args.copy()
    args[1] = str(source)
    assert main(args + ['--overwrite']) == 1
    assert source.read_bytes() == original
    assert not (sample.output / 'result.run.json').exists()


@pytest.mark.parametrize('contents', ['sample_id\tfastq\tfastq\na\tb\tc\n', 'sample_id\tfastq\na\n', 'sample_id\tfastq\na\tb\na\tc\n'])
def test_malformed_batch_is_rejected_before_running(tmp_path, contents):
    manifest = tmp_path / 'cells.tsv'
    manifest.write_text(contents)
    with pytest.raises(ValueError):
        read_batch(manifest, parse_args([]))


def test_publish_failure_rolls_back_all_previous_files(tmp_path, monkeypatch):
    stage, dest = tmp_path / 'stage', tmp_path / 'dest'
    stage.mkdir()
    dest.mkdir()
    for name in ('one', 'two'):
        (stage / name).write_text('new')
        (dest / name).write_text('old')
    rename = pipeline.os.replace
    def fail_second(source, target):
        if Path(source) == stage / 'two':
            raise OSError('simulated disk error')
        return rename(source, target)
    monkeypatch.setattr(pipeline.os, 'replace', fail_second)
    with pytest.raises(OSError):
        pipeline.publish_bundle(stage, dest, ['one', 'two'])
    assert (dest / 'one').read_text() == (dest / 'two').read_text() == 'old'


def test_mouse_tcr_chain_order_is_derived_from_references():
    db = Path(__file__).resolve().parents[1] / 'db' / 'mouse' / 'TCR'
    mappings = {name: SimpleNamespace(candidates=[]) for name in ('hv', 'hc', 'lv', 'lc')}
    assert pipeline._tcr_chain_names({'db_path': str(db)}, mappings) == ('TCR_1', 'TCR_2')
    for region, locus in [('hv', 'TRB'), ('hc', 'TRB'), ('lv', 'TRA'), ('lc', 'TRA')]:
        mappings[region].candidates = [{'sequence': 'ACGT', 'count': 2, 'reference': locus+'V1', 'locus': locus, 'mapq': 42}]
    assert pipeline._tcr_chain_names({'db_path': str(db)}, mappings) == ('TRB', 'TRA')


def test_failed_rollback_retains_previous_files_outside_staging(tmp_path, monkeypatch):
    stage, dest = tmp_path / 'stage', tmp_path / 'dest'
    stage.mkdir()
    dest.mkdir()
    for name in ('one', 'two'):
        (stage / name).write_text('new')
        (dest / name).write_text('old')
    rename = pipeline.os.replace
    def fail_install_and_restore(source, target):
        source = Path(source)
        if source == stage / 'two' or source.parent.name.startswith('.basic-recovery-'):
            raise OSError('persistent write error')
        return rename(source, target)
    monkeypatch.setattr(pipeline.os, 'replace', fail_install_and_restore)
    with pytest.raises(RuntimeError, match='recoverable'):
        pipeline.publish_bundle(stage, dest, ['one', 'two'])
    backup = next(dest.glob('.basic-recovery-*'))
    assert (backup / 'one').read_text() == (backup / 'two').read_text() == 'old'


def test_incompatible_loci_are_not_stitched_together():
    def candidate(locus):
        return {'sequence': 'ACGT', 'locus': locus, 'reference': locus + '1', 'mapq': 42, 'count': 3}
    pairs = pipeline._seed_pairs([candidate('IGK')], [candidate('IGL')], 1)
    assert all(v is None or c is None for v, c in pairs)


def test_seed_pairs_preserve_all_four_tcr_loci_with_one_anchor_per_locus():
    def candidates():
        return [{'sequence': sequence, 'locus': locus, 'reference': locus + '1', 'mapq': 42, 'count': 3}
                for locus, sequence in zip(('TRA', 'TRB', 'TRG', 'TRD'), ('AAAA', 'CCCC', 'GGGG', 'TTTT'))]
    pairs = pipeline._seed_pairs(candidates(), candidates(), 1)
    assert len(pairs) == 4
    assert {v['locus'] for v, c in pairs} == {'TRA', 'TRB', 'TRG', 'TRD'}
    assert all(v['locus'] == c['locus'] for v, c in pairs)


def test_chain_summary_locus_follows_best_retained_assembly(tmp_path):
    args = parse_args([])
    args.update(type='TCR', min_overlap=2, db_path=str(tmp_path))
    def candidate(locus, sequence, count):
        return dict(locus=locus, sequence=sequence, count=count, reference=locus + '1', mapq=42)
    mappings = {
        'hv': SimpleNamespace(candidates=[candidate('TRB', 'AAAA', 10), candidate('TRD', 'ACGT', 3)]),
        'hc': SimpleNamespace(candidates=[candidate('TRB', 'CCCC', 10), candidate('TRD', 'ACGT', 3)]),
        'lv': SimpleNamespace(candidates=[candidate('TRA', 'GTGT', 3)]),
        'lc': SimpleNamespace(candidates=[candidate('TRA', 'ACAC', 3)]),
    }
    pool = SimpleNamespace(counts={'AAAA': 10, 'CCCC': 10, 'ACGT': 3, 'GTGT': 3, 'ACAC': 3},
                           quality_weights={}, mate_links={})
    chains, records = pipeline._assemble_results(args, pool, mappings)
    # Higher-support TRB seeds cannot connect. Only the TRD seed pair completes.
    assert chains[0]['chain'] == 'TRD'
    assert chains[0]['status'] == 'complete'
    assert chains[0]['region_group'] == 'hv'
    assert {record['locus'] for record in records if record['basic_region_group'] == 'hv'} == {'TRD'}
    summary = tmp_path / 'summary.tsv'
    pipeline._summary(summary, chains, records)
    with summary.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert rows[0]['chain'] == 'TRD' and rows[0]['status'] == 'complete'


def test_allele_mapping_does_not_split_identical_seed_support():
    candidates = [{'sequence': 'ACGT', 'locus': 'IGH', 'reference': r, 'mapq': 0, 'count': 2}
                  for r in ('IGHV1', 'IGHV2')]
    aggregate = pipeline._aggregate(candidates)
    assert len(aggregate) == 1 and aggregate[0]['count'] == 4


def test_cli_check_does_not_create_outputs(sample):
    assert main(sample.args + ['--check']) == 0
    assert not sample.output.exists()


def test_keep_temp_records_prepared_reads_and_logs(sample):
    assert main(sample.args + ['--keep-temp']) == 0
    manifest = json.loads((sample.output / 'result.run.json').read_text())
    work = Path(manifest['temporary_directory'])
    assert (work / 'reads.fastq').is_file()
    assert len(list(work.glob('*.stderr.log'))) == 4


def test_assembly_exception_fails_explicitly_without_worker_wait(sample, monkeypatch):
    def crash(*args, **kwargs):
        raise ValueError('synthetic assembly failure')
    monkeypatch.setattr(pipeline, 'assemble_chain', crash)
    assert main(sample.args) == 1
    failure = json.loads((sample.output / 'result.failed.json').read_text())
    assert failure['error_type'] == 'ValueError'
    assert not (sample.output / 'result.fasta').exists()


def test_strict_mode_does_not_publish_partial_alternatives(sample, monkeypatch):
    original = pipeline._seed_pairs
    def pairs(*args, **kwargs):
        first = original(*args, **kwargs)[0]
        alternative = (dict(first[0], sequence='ACGTACGTACGTACGT', count=1), first[1])
        return [first, alternative]
    def assemble(start, *args, **kwargs):
        partial = start == 'ACGTACGTACGTACGT'
        return {'status': 'partial' if partial else 'complete', 'contigs': [
            {'sequence': start, 'region': 'variable' if partial else 'full', 'ambiguous': False}], 'diagnostics': {}}
    monkeypatch.setattr(pipeline, '_seed_pairs', pairs)
    monkeypatch.setattr(pipeline, 'assemble_chain', assemble)
    assert main(sample.args + ['--max-candidates', '2']) == 0
    records = read_airr(sample.output / 'result.airr.tsv')
    assert len(records) == 2 and all(r['basic_status'] == 'complete' for r in records)


def test_unmatched_locus_is_retained_as_partial_candidate():
    def candidate(locus):
        return {'sequence': 'ACGT', 'locus': locus, 'reference': locus+'1', 'mapq': 42, 'count': 3}
    pairs = pipeline._seed_pairs([candidate('IGK'), candidate('IGL')], [candidate('IGK')], 1)
    assert any(v and v['locus'] == 'IGL' and c is None for v, c in pairs)


def test_cli_rejects_mixed_input_modes(sample):
    assert main(sample.args + ['-PE_1', str(sample.reads), '-PE_2', str(sample.reads)]) == 1


def test_related_species_keeps_target_and_seed_species_separate(sample):
    assert main(sample.args + ['-g', 'nonmodel', '--seed-species', 'human', '--local-seeding']) == 0
    manifest = json.loads((sample.output / 'result.run.json').read_text())
    assert manifest['target_species'] == 'nonmodel'
    assert manifest['seed_species'] == 'human'
    assert manifest['experimental'] is True
    assert '--very-sensitive-local' in manifest['mapping']['hv']['command']
    assert all(r['basic_reference_transfer'] == 'T' for r in read_airr(sample.output / 'result.airr.tsv'))


def test_annotation_preflight_does_not_claim_missing_executable_valid(sample):
    assert main(sample.args + ['--check', '--igblast', str(sample.tmp / 'absent'), '--germline-v', 'V', '--germline-j', 'J']) == 1
    assert not sample.output.exists()


def test_output_cannot_replace_raw_input_even_with_overwrite(sample):
    original = sample.reads.read_bytes()
    target = sample.tmp / 'result.fasta'
    sample.reads.rename(target)
    args = sample.args.copy()
    args[1] = str(target)
    args[args.index('-o') + 1] = str(sample.tmp)
    assert main(args + ['--overwrite']) == 1
    assert target.read_bytes() == original


def test_batch_sample_ids_cannot_alias_by_case(tmp_path):
    manifest = tmp_path / 'cells.tsv'
    manifest.write_text('sample_id\tfastq\nCellA\tone.fastq\ncella\ttwo.fastq\n')
    with pytest.raises(ValueError, match='duplicate'):
        read_batch(manifest, parse_args([]))


def test_cli_module_and_legacy_entrypoints_agree():
    root = Path(__file__).resolve().parents[1]
    for command in ([sys.executable, str(root / 'BASIC.py')], [sys.executable, '-m', 'basic_assembler']):
        result = subprocess.run(command + ['--version'], cwd=root, text=True, capture_output=True)
        from basic_assembler import __version__
        assert result.returncode == 0 and __version__ in result.stdout
