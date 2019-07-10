import BASIC
import pytest
import os
from subprocess import CalledProcessError
from . import utils


def test_complement():
    assert BASIC.complement('A') == 'T'
    assert BASIC.complement('T') == 'A'
    assert BASIC.complement('G') == 'C'
    assert BASIC.complement('C') == 'G'


def test_reverse_complement():
    assert BASIC.revcom('ATGC') == 'GCAT'


class TestConfig(object):

    def setup_method(self):
        self.repo_path = os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)
            ))

        self.args = {
            'type': 'BCR',
            'num_threads': 2,
            'name': 'result',
            'FASTQ': '',
            'LEFT': os.path.join(
                self.repo_path,
                'examples/H8_AW1/R1_100K_reads.fastq.gz'),
            'RIGHT': os.path.join(
                self.repo_path,
                'examples/H8_AW1/R2_100K_reads.fastq.gz'),
            'genome': 'human',
            'bowtie': utils.find_bowtie2(),
            'tmpdir': './',
            'output_location': './',
            'allow_partial': False,
            'VERBOSE': False,
            'database_path': '',
        }

    def test_default_args_with_ex_fastqs(self):
        args_out = BASIC.config_args(self.args)
        assert args_out['read_type'] == 'paired'

    def test_single_end_fastq(self):
        self.args['LEFT'] = ''
        self.args['RIGHT'] = ''
        self.args['FASTQ'] = 'examples/cell_A1/ERR1421621_1M_reads.fastq.gz'
        args_out = BASIC.config_args(self.args)
        assert args_out['read_type'] == 'single'

    def test_missing_seqs(self):
        self.args['LEFT'] = ''
        self.args['RIGHT'] = ''
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_missing_left_of_paired_seqs(self):
        self.args['LEFT'] = ''
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_missing_right_of_paired_seqs(self):
        self.args['RIGHT'] = ''
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_invalid_assembly_type(self):
        self.args['type'] = 'BCX'
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_invalid_genome(self):
        self.args['genome'] = 'HUMAX'
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_invalid_db_path(self):
        self.args['database_path'] = os.path.abspath(__file__)
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_invalid_output_location(self):
        self.args['output_location'] = '/nopermission'
        with pytest.raises(SystemExit):
            BASIC.config_args(self.args)

    def test_wrong_bowtie2_path(self):
        self.args['bowtie'] = 'nonexistent'
        with pytest.raises(EnvironmentError):
            BASIC.config_args(self.args)

    def test_bowtie2_passing_directory(self):
        self.args['bowtie'] = os.path.dirname(self.args['bowtie'])
        args_out = BASIC.config_args(self.args)
        assert args_out['read_type'] == 'paired'
