# -*- coding: utf-8 -*-

from typing import NewType, Tuple, Sequence, ValuesView

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

COORDS = 0
OVL_LEN = 1


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

    def __delitem__(self, key):
        del self._matrix[key]
    # end def __delitem__

    def get_overlaps_by_index(self, index: ContigIndex) -> Sequence[Overlap]:
        return tuple(
            filter(
                lambda x: x[COORDS][0] == index,
                self._matrix.items()
            )
        )
    # end def get_overlaps_by_key

    def get_keys_by_index(self, index: ContigIndex) -> Sequence[Overlap]:
        return tuple(
            filter(
                lambda x: x[0] == index,
                self._matrix.keys()
            )
        )
    # end def get_overlaps_by_key

    def items(self) -> ValuesView:
        return self._matrix.items()
    # end def values
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
    SmE, EmS, SmS, EmE = matrices

    # For comments for loop below, let x denote to start of a contig and y -- to end of a contig.
    # On the second iteration it will be vice versa: x will be end of a contig and y -- start.
    # But for clarity, in comments we will write explicitly: "start" for x and "end" for y.

    XmY: MatchMatrix
    YmX: MatchMatrix
    XmX: MatchMatrix
    for XmY, YmX, XmX in zip((SmE, EmS), (EmS, SmE), (SmS, EmE)):
        i: ContigIndex
        for i in range(len(contig_collection)):

            max_overlap = _get_largest_ovl(i, XmY, XmX)

            if max_overlap[OVL_LEN] != 0:
                # Find number of contig, to which current start is adjasent
                mate_idx = max_overlap[COORDS][1]
                # Amend genome length -- substract length of the overlap
                #   multiplied by minimum of multiplicity of adjacent contigs
                expexted_genome_size -= max_overlap[OVL_LEN]\
                                        * min(
                                              contig_collection[i].multplty,
                                              contig_collection[mate_idx].multplty
                                          )

                _clear_matrices(i, mate_idx, XmY, YmX, XmX)
            # end if
        # end for
    # end for

    return int(expexted_genome_size)

# end def _trace_back


def _get_largest_ovl(i, XmY, XmX):
    XmY_ovls = XmY.get_overlaps_by_index(i)
    XmX_ovls = XmX.get_overlaps_by_index(i)
    if len(XmY_ovls) != 0 or len(XmX_ovls) != 0:
        max_overlap = max(XmY_ovls + XmX_ovls,
                          key=lambda x: x[1]
                      )
    else:
        max_overlap = ((None, None), 0)
    # end if
    return max_overlap
# end def _get_largest_ovl


def _clear_matrices(i, mate_idx, XmY, YmX, XmX):

    # Set all elements in XmY, YmX and XmX, which overlaps with current start.
    # I.e. not only maximum one -- all.
    # XmY[j][i] has nothing to do with current start -- it should be left >0, if it is.
    # That is why we actually need XmX matrices.
    for coords in XmY.get_keys_by_index(i):
        del XmY[coords]
    # end for

    for matrix in (YmX, XmX):
        for coords in matrix.get_keys_by_index(i):
            del matrix[coords]
        # end for
        for coords in matrix.get_keys_by_index(mate_idx):
            del matrix[coords]
        # end for
    # end for
# end def _clear_matrices
