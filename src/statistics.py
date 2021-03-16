# -*- coding: utf-8 -*-

import typing

from src.contigs import Contig
from src.overlaps import START, RCSTART, END, RCEND, Overlap, OverlapCollection


ContigCollection = typing.MutableSequence[Contig]


_START2START_MATCH = {START, RCSTART}
_END2END_MATCH = {END, RCEND}


def calc_sum_contig_lengths(contig_collection: ContigCollection):
    return sum(map(lambda x: x.length, contig_collection))
# end def calc_sum_contig_lengths


def get_min_coverage(contig_collection: ContigCollection):
    return min(map(lambda x: x.cov, contig_collection))
# end def get_min_coverage


def get_max_coverage(contig_collection: ContigCollection):
    return max(map(lambda x: x.cov, contig_collection))
# end def get_min_coverage


def calc_mean_coverage(contig_collection: ContigCollection):
    mean_cov: float = sum(map(lambda x: x.cov, contig_collection))\
                      / len(contig_collection)
    return round(mean_cov, 2)
# end def get_min_coverage


def calc_lq_coef(contig_collection: ContigCollection, overlap_collection: OverlapCollection):
    total_dead_ends: int = 0
    total_termini: int = 0

    i: int
    contig: Contig
    num_termini: int
    num_start_matches: int
    num_end_matches: int
    num_dead_ends: int

    for i, contig in enumerate(contig_collection):
        num_termini = int(2 * contig.multplty)
        total_termini += num_termini

        num_start_matches = len(tuple(filter(_is_start_match, overlap_collection[i])))
        num_end_matches   = len(tuple(filter(_is_end_match, overlap_collection[i])))

        num_dead_ends = int(2 * contig.multplty)\
                        - min(contig.multplty, num_start_matches)\
                        - min(contig.multplty, num_end_matches)

        total_dead_ends += num_dead_ends
    # end for

    lq_coef: int = (1 - total_dead_ends / total_termini) * 100

    return round(lq_coef, 2)
# end def calc_lq_coef


def _is_start_match(ovl: Overlap):
    return (ovl.terminus1 == START and ovl.terminus2 == END)\
        or {ovl.terminus1, ovl.terminus2} == _START2START_MATCH
# end def def _is_start_match


def _is_end_match(ovl: Overlap):
    return (ovl.terminus1 == END and ovl.terminus2 == START)\
        or {ovl.terminus1, ovl.terminus2} == _END2END_MATCH
# end def def _is_end_match
