# -*- coding: utf-8 -*-

def find_overlap_s2s(seq1: str, seq2: str, mink: int, maxk: int) -> int:
    # Function searches for identity between starts of seq1 and seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;
    # :param seq2: another sequence;
    # :param mink: minimun overlap;
    # :param maxk: maximum overlap;
    #
    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    i: int = mink
    overlap: int = 0

    while seq1[:i] == seq2[:i] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap_s2s


def find_overlap_e2s(seq1: str, seq2: str, mink: int, maxk: int) -> int:
    # Function searches for identity between end of seq1 and start of seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;
    # :param seq2: another sequence;
    # :param mink: minimun overlap;
    # :param maxk: maximum overlap;
    #
    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    ovl_len: int = 0

    i: int
    for i in range(mink, maxk+1):
        if seq1[-i:] == seq2[:i]:
            ovl_len = i
        # end if
    # end for

    return ovl_len
# end def find_overlap_e2s


def find_overlap_e2e(seq1: str, seq2: str, mink: int, maxk: int) -> int:
    # Function searches for identity between ends of seq1 and seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;
    # :param seq2: another sequence;
    # :param mink: minimun overlap;
    # :param maxk: maximum overlap;
    #
    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    i: int = mink
    overlap: int = 0

    while seq1[-i:] == seq2[-i:] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap_e2e
