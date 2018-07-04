import BASIC


def test_complement():
    assert BASIC.complement('A') == 'T'
    assert BASIC.complement('T') == 'A'
    assert BASIC.complement('G') == 'C'
    assert BASIC.complement('C') == 'G'


def test_reverse_complement():
    assert BASIC.revcom('ATGC') == 'GCAT'
