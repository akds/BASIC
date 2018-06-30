from Bio import SeqIO


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


def test_se_human_bcr():
    expected = read_fasta('examples/cell_A1/result.fasta')
    actual = read_fasta('SE_test.fasta')

    for chain_id, seq in expected.items():
        assert actual[chain_id] == seq
