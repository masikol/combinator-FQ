# -*- coding: utf-8 -*-

from typing import NewType, Tuple

from src.contigs import ContigIndex, ContigCollection
from src.overlaps import Overlap, OverlapCollection
from src.overlaps import Terminus, START, RCSTART, END, RCEND
from src.combinator_statistics import calc_sum_contig_lengths

# Matrices for amending expected genome length:
#
# Explanation:
# For example, we have matrix SmE (maxtrix of starts matching ends):
#
# 0 0 0
# 5 0 0
# 0 0 7
#
# It is matrix 3x3, therefore, we have 3 contigs.
# SmE[2][2] == 7, therefore contig number 2 (0-based index) is circular with overlap of 7 b.p.
# SmE[0][1] == 5, therefore start of contig #0 matches end of contig #1.


OverlapCoords = NewType('OverlapCoords', Tuple[int, int])


class MatchMatrix:

    def __init__(self):
        self._matrix = dict()
    # end def __init__

    def __getitem__(self, coords: OverlapCoords) -> int:
        try:
            item_to_return = self._matrix[coords]
        except KeyError:
            item_to_return = 0
        # end try
        return item_to_return
    # end def __getitem__

    def __setitem__(self, coords: OverlapCoords, ovl_len: int) -> None:
        self._matrix[coords] = ovl_len
    # end def __setitem__

    def __repr__(self) -> str:
        return str(self._matrix)
    # end def __repr__
# end class MatchMatrix


def calc_exp_genome_size(contig_collection: ContigCollection,
                         overlap_collection: OverlapCollection) -> int:
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;

    matrices = _fill_matrices(overlap_collection)
    expeted_genome_size = _trace_back(matrices, contig_collection)

    return expeted_genome_size
# end def calc_exp_genome_size


def _fill_matrices(overlap_collection: OverlapCollection) -> Tuple[MatchMatrix, MatchMatrix,
                                                                   MatchMatrix, MatchMatrix]:
    # Matrix for starts matching ends:
    SmE: MatchMatrix = MatchMatrix()
    # Matrix for ends matching starts:
    EmS: MatchMatrix = MatchMatrix()
    # Matrix for starts matching starts (reverse-complement matching):
    SmS: MatchMatrix = MatchMatrix()
    # Matrix for ends matching ends (reverse-complement matching):
    EmE: MatchMatrix = MatchMatrix()

    s2e_match: Tuple[Terminus, Terminus] = (START, END)
    e2s_match: Tuple[Terminus, Terminus] = (END, START)
    s2s_match: Tuple[Terminus, Terminus] = (START, RCSTART)
    e2e_match: Tuple[Terminus, Terminus] = (END, RCEND)

    key: ContigIndex
    for key in range(len(overlap_collection)):
        overlaps = overlap_collection[key]

        if len(overlaps) != 0:

            ovl: Overlap
            for ovl in overlaps:
                termini: Tuple[Terminus, Terminus] = (ovl.terminus_i, ovl.terminus_j)

                if termini == s2e_match:
                    SmE[(ovl.contig_i, ovl.contig_j)] = ovl.ovl_len

                elif termini == e2s_match:
                    EmS[(ovl.contig_i, ovl.contig_j)] = ovl.ovl_len

                elif termini == s2s_match:
                    SmS[(ovl.contig_i, ovl.contig_j)] = ovl.ovl_len

                elif termini == e2e_match:
                    EmE[(ovl.contig_i, ovl.contig_j)] = ovl.ovl_len
                # end if
            # end for
        # end if
    # end for

    return SmE, EmS, SmS, EmE
# end def _fill_matrices

def _trace_back(
    matrices: Tuple[MatchMatrix, MatchMatrix, MatchMatrix, MatchMatrix],
    contig_collection: ContigCollection
) -> int:

    expexted_genome_size: int = calc_sum_contig_lengths(contig_collection)


# end def _trace_back
