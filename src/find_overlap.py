# -*- coding: utf-8 -*-

def find_overlap_s2s(seq1, seq2, mink, maxk):
    # Function searches for identity between starts of seq1 and seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;      :type seq1: str;
    # :param seq2: another sequence;  :type seq2: str;
    # :param mink: minimun overlap;   :type mink: int;
    # :param maxk: maximum overlap;   :type maxk: int;

    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    i = mink
    overlap = 0

    while seq1[:i] == seq2[:i] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap_s2s


def find_overlap_e2s(seq1, seq2, mink, maxk):
    # Function searches for identity between end of seq1 and start of seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;      :type seq1: str;
    # :param seq2: another sequence;  :type seq2: str;
    # :param mink: minimun overlap;   :type mink: int;
    # :param maxk: maximum overlap;   :type maxk: int;
    #
    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    for i in range(mink, maxk+1):
        if seq1[-i:] == seq2[:i]:
            return i
        # end if
    # end for

    return 0
# end def find_overlap_e2s


def find_overlap_e2e(seq1, seq2, mink, maxk):
    # Function searches for identity between ends of seq1 and seq2.
    # Function regards overlap of length [mink, maxk].
    #
    # :param seq1: one sequence;      :type seq1: str;
    # :param seq2: another sequence;  :type seq2: str;
    # :param mink: minimun overlap;   :type mink: int;
    # :param maxk: maximum overlap;   :type maxk: int;
    #
    # Returns 0 if overlap is less than 'mink' and
    #   length of overlap (which is <= maxk) otherwise.

    i = mink
    overlap = 0

    while seq1[-i:] == seq2[-i:] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap_e2e
