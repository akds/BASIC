from Bio import SeqIO
import os
import subprocess
import pytest


def find_bowtie2_dir():
    """ Looks for bowtie2 in PATH if specific environmental variable does
        not exist
    """

    try:
        return os.environ['bowtie2']
    except KeyError:
        stdout, stderr = subprocess.Popen(["which", "bowtie2"],
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE).communicate()
        if stdout:
            return os.path.dirname(stdout.strip())
        else:
            print('Bowtie2 not found. Please add it to $PATH or export the '
                  'environmental variable "bowtie2" specifying its directory')
            exit(1)


def read_fasta(infile):
    """ Returns fasta as a dict of {chain: sequence}
    """
    recs = {}
    for record in SeqIO.parse(infile, 'fasta'):
        # extract chain from header e.g. >cell_id=result;heavy_chain;BASIC
        chain = record.id.split(';')[1]
        seq = str(record.seq)
        recs[chain] = seq
    return recs


def run_basic(cmd):
    bowtie2_path = find_bowtie2_dir()
    cmd += " -b {}".format(bowtie2_path)

    subprocess.check_call(cmd.split(' '))

def test_se_human_bcr():

    cmd = ("python BASIC.py -SE examples/cell_A1/ERR1421621_1M_reads.fastq.gz "
           "-g human -i BCR -n SE_test")

    run_basic(cmd)

    expected = read_fasta('examples/cell_A1/result.fasta')
    actual = read_fasta('SE_test.fasta')

    for chain_id, seq in expected.items():
        assert actual[chain_id] == seq


def test_pe_human_bcr():

    # run BASIC
    cmd = ("python BASIC.py "
           "-PE_1 examples/H8_AW1/R1_100K_reads.fastq.gz "
           "-PE_2 examples/H8_AW1/R2_100K_reads.fastq.gz "
           "-g human -i BCR -n PE_test")

    run_basic(cmd)

    sanger = read_fasta('examples/H8_AW1/sanger.fasta')
    actual = read_fasta('PE_test.fasta')

    # check assembly is inclusive of Sanger result (assembly will contain
    # additional sequence e.g. V leader)
    for chain_id, seq in sanger.items():
        assert seq in actual[chain_id]


def test_empty_fail_without_partial_arg():

    subprocess.check_call(["touch", "empty_file.fasta"])

    cmd = ("python BASIC.py -SE empty_file.fasta "
           "-g human -i BCR -n no_partial")

    with pytest.raises(subprocess.CalledProcessError):
        run_basic(cmd)


def test_empty_success_on_partial_arg():

    # an empty fasta as input should pass with the partial argument
    subprocess.check_call(["touch", "empty_file.fasta"])

    cmd = ("python BASIC.py -SE empty_file.fasta "
           "-g human -i BCR -n partial_success -a")

    run_basic(cmd)

    # empty result file
    assert os.stat('partial_success.fasta').st_size == 0
