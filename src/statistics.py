# -*- coding: utf-8 -*-

from typing import Collection

from src.contigs import Contig, ContigCollection
from src.overlaps import Terminus, START, RCSTART, END, RCEND
from src.overlaps import Overlap, OverlapCollection


class CoverageCalculator:

    def __init__(self, contig_collection: ContigCollection):
        self.coverages: Collection = self._filter_non_none_covs(contig_collection)
        self.num_contigs: int = len(contig_collection)
    # end def __init__

    def get_min_coverage(self) -> float:

        min_cov: float
        if len(self.coverages) == 0:
            min_cov = None
        else:
            min_cov = min(self.coverages)
        # end if

        return min_cov
    # end def get_min_coverage


    def get_max_coverage(self) -> float:

        max_cov: float
        if len(self.coverages) == 0:
            max_cov = None
        else:
            max_cov = max(self.coverages)
        # end if

        return max_cov
    # end def get_min_coverage


    def calc_mean_coverage(self) -> float:

        mean_cov: float
        if len(self.coverages) == 0:
            max_cov = None
        else:
            mean_cov = sum(self.coverages) / self.num_contigs
            mean_cov = round(mean_cov, 2)
        # end if

        return mean_cov
    # end def get_min_coverage


    def _filter_non_none_covs(self,
                              contig_collection: ContigCollection) -> Collection[float]:
        return tuple(
            filter(
                lambda x: not x is None,
                map(
                    lambda x: x.cov, contig_collection
                )
            )
        )
    # end def _filter_non_none_covs
# end class CoverageCalculator


_START2START_MATCH: set = {START, RCSTART}
_END2END_MATCH:     set = {END, RCEND}


def calc_sum_contig_lengths(contig_collection: ContigCollection) -> int:
    return sum(
        map(lambda x: x.length, contig_collection)
    )
# end def calc_sum_contig_lengths


def calc_lq_coef(contig_collection: ContigCollection,
                 overlap_collection: OverlapCollection) -> float:
    total_termini: int = 0
    total_dead_ends: int = 0

    i: int
    contig: Contig

    for i, contig in enumerate(contig_collection):

        round_multplty: int = round(contig.multplty)
        num_termini: int = int(2 * round_multplty)

        num_start_matches: int = len(tuple(
            filter(is_start_match, overlap_collection[i])
        ))
        num_end_matches:   int = len(tuple(
            filter(is_end_match, overlap_collection[i])
        ))

        num_dead_ends: int = int(2 * round_multplty)\
                             - min(round_multplty, num_start_matches)\
                             - min(round_multplty, num_end_matches)

        total_termini += num_termini
        total_dead_ends += num_dead_ends
    # end for

    lq_coef: float = (1 - total_dead_ends / total_termini) * 100.0

    return round(lq_coef, 2)
# end def calc_lq_coef


def is_start_match(ovl: Overlap) -> bool:
    return (ovl.terminus1 == START and ovl.terminus2 == END)\
           or {ovl.terminus1, ovl.terminus2} == _START2START_MATCH
# end def def is_start_match


def is_end_match(ovl: Overlap) -> bool:
    return (ovl.terminus1 == END and ovl.terminus2 == START)\
           or {ovl.terminus1, ovl.terminus2} == _END2END_MATCH
# end def def is_end_match
