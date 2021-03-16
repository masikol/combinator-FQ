# -*- coding: utf-8 -*-

from src.overlaps import START, RCSTART, END, RCEND


_START2START_MATCH = {START, RCSTART}
_END2END_MATCH = {END, RCEND}
_START2END_MATCH = {END, START}


def calc_sum_contig_lengths(contig_collection):
    return sum(map(lambda x: x.length, contig_collection))
# end def calc_sum_contig_lengths


def get_min_coverage(contig_collection):
    return min(map(lambda x: x.cov, contig_collection))
# end def get_min_coverage


def get_max_coverage(contig_collection):
    return max(map(lambda x: x.cov, contig_collection))
# end def get_min_coverage


def calc_mean_coverage(contig_collection):
    return round(
        sum(map(lambda x: x.cov, contig_collection))
            / len(contig_collection),
        2
    )
# end def get_min_coverage


def calc_LQ_coef(contig_collection, overlap_collection):
    total_dead_ends = 0
    total_termini = 0

    for i, contig in enumerate(contig_collection):
        num_termini = int(2 * contig.multplty)
        total_termini += num_termini

        start_matches_num = len(tuple(filter(_is_start_match, overlap_collection[i])))
        end_matches_num   = len(tuple(filter(_is_end_match, overlap_collection[i])))
        num_dead_ends = int(2 * contig.multplty)\
                        - min(contig.multplty, start_matches_num)\
                        - min(contig.multplty, end_matches_num)
        total_dead_ends += num_dead_ends
    # end for

    lq = (1 - total_dead_ends / total_termini) * 100

    return round(lq, 2)
# end def calc_LQ_coef


def _is_start_match(ovl):
    return (ovl.terminus1 == START and ovl.terminus2 == END)\
        or {ovl.terminus1, ovl.terminus2} == _START2START_MATCH
# end def def _is_start_match


def _is_end_match(ovl):
    return (ovl.terminus1 == END and ovl.terminus2 == START)\
        or {ovl.terminus1, ovl.terminus2} == _END2END_MATCH
# end def def _is_end_match
