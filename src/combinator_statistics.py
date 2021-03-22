# -*- coding: utf-8 -*-

from typing import Collection
from statistics import mean, median

from src.contigs import Contig, ContigCollection
from src.overlaps import START, RCSTART, END, RCEND
from src.overlaps import Overlap, OverlapCollection


class CoverageCalculator:
    # Class for calcuating statistics associated with coverage.

    def __init__(self, contig_collection: ContigCollection):
        # :param contig_collection: instance of ContigCollection returned by
        #   `src.contigs.get_contig_collection` function;
        self.coverages: Collection = self._filter_non_none_covs(contig_collection)
    # end def __init__

    def get_min_coverage(self) -> float:
        # Method for obtaining minimum coverage value.

        min_cov: float
        if len(self.coverages) == 0:
            min_cov = None # no coverage values are available
        else:
            min_cov = min(self.coverages)
        # end if

        return min_cov
    # end def get_min_coverage


    def get_max_coverage(self) -> float:
        # Method for obtaining maximum coverage value.

        max_cov: float
        if len(self.coverages) == 0:
            max_cov = None # no coverage values are available
        else:
            max_cov = max(self.coverages)
        # end if

        return max_cov
    # end def get_min_coverage


    def calc_mean_coverage(self) -> float:
        # Method for calculating mean coverage.

        mean_cov: float
        if len(self.coverages) == 0:
            mean_cov = None # no coverage values are available
        else:
            mean_cov = mean(self.coverages)
            mean_cov = round(mean_cov, 2)
        # end if

        return mean_cov
    # end def get_min_coverage

    def calc_median_coverage(self) -> float:
        # Method for calculating median coverage.

        median_cov: float
        if len(self.coverages) == 0:
            median_cov = None # no coverage values are available
        else:
            median_cov = median(self.coverages)
            median_cov = round(median_cov, 2)
        # end if

        return median_cov
    # end def get_min_coverage

    @staticmethod
    def _filter_non_none_covs(contig_collection: ContigCollection) -> Collection[float]:
        # Method selects contigs with non-None coverage value from a `contig_collection`.
        #
        # :param contig_collection: instance of ContigCollection returned by
        #   `src.contigs.get_contig_collection` function;
        # Returns tuple of coverage values.

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


# Only start-rcstart matches is adjacency-associated fro start-to-start matches
_START2START_MATCH: set = {START, RCSTART}
# Only end-rcend matches is adjacency-associated fro end-to-end matches
_END2END_MATCH:     set = {END, RCEND}


def calc_sum_contig_lengths(contig_collection: ContigCollection) -> int:
    # Function for summarizing contigs lengths.
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    return sum(
        map(lambda x: x.length, contig_collection)
    )
# end def calc_sum_contig_lengths


def calc_lq_coef(contig_collection: ContigCollection,
                 overlap_collection: OverlapCollection) -> float:
    # Function calculates LQ-coefficient for given contig set.
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;

    # Total number of termini taking account of multiplicity
    total_termini: int = 0
    # Total number of dead ends taking account of multiplicity
    total_dead_ends: int = 0

    i: int
    contig: Contig
    # Iterate over contigs
    for i, contig in enumerate(contig_collection):

        # Get rounded multiplicity
        round_multplty: int = round(contig.multplty)
        # Calculate number of termini taking account of multiplicity
        num_termini: int = int(2 * round_multplty)

        # Cound overlaps associated with start
        num_start_matches: int = len(tuple(
            filter(is_start_match, overlap_collection[i])
        ))
        # Cound overlaps associated with end
        num_end_matches:   int = len(tuple(
            filter(is_end_match, overlap_collection[i])
        ))

        # Calculate number of dead ends of the current contig
        num_dead_ends: int = int(2 * round_multplty)\
                             - min(round_multplty, num_start_matches)\
                             - min(round_multplty, num_end_matches)

        # Increment `total_termini` and `total_dead_ends`
        total_termini += num_termini
        total_dead_ends += num_dead_ends
    # end for

    # Calculate the LQ coefficient
    lq_coef: float = (1 - total_dead_ends / total_termini) * 100.0

    return round(lq_coef, 2)
# end def calc_lq_coef


def is_start_match(ovl: Overlap) -> bool:
    # Function returns True if overlap `ovl` is associated with start.
    return (ovl.terminus1 == START and ovl.terminus2 == END)\
           or {ovl.terminus1, ovl.terminus2} == _START2START_MATCH
# end def def is_start_match


def is_end_match(ovl: Overlap) -> bool:
    # Function returns True if overlap `ovl` is associated with end.
    return (ovl.terminus1 == END and ovl.terminus2 == START)\
           or {ovl.terminus1, ovl.terminus2} == _END2END_MATCH
# end def def is_end_match
