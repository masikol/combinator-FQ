# -*- encoding: utf-8 -*-

from typing import Collection, Sequence, Callable
from statistics import mean, median

from src.contigs import ContigIndex, Contig, ContigCollection
from src.overlaps import START, RCSTART, END, RCEND
from src.overlaps import Overlap, OverlapCollection


class CoverageCalculator:
    # Class for calcuating statistics associated with coverage.

    def __init__(self, contig_collection: ContigCollection):
        # :param contig_collection: instance of ContigCollection returned by
        #   `src.contigs.get_contig_collection` function;
        self._coverages: Collection[float] = self._filter_non_none_covs(contig_collection)
    # end def __init__

    def get_min_coverage(self) -> float:
        # Method for obtaining minimum coverage value.

        min_cov: float
        if len(self._coverages) == 0:
            min_cov = None # no coverage values are available
        else:
            min_cov = min(self._coverages)
        # end if

        return min_cov
    # end def get_min_coverage


    def get_max_coverage(self) -> float:
        # Method for obtaining maximum coverage value.

        max_cov: float
        if len(self._coverages) == 0:
            max_cov = None # no coverage values are available
        else:
            max_cov = max(self._coverages)
        # end if

        return max_cov
    # end def get_min_coverage


    def calc_mean_coverage(self) -> float:
        # Method for calculating mean coverage.

        mean_cov: float
        if len(self._coverages) == 0:
            mean_cov = None # no coverage values are available
        else:
            mean_cov = mean(self._coverages)
            mean_cov = round(mean_cov, 2)
        # end if

        return mean_cov
    # end def get_min_coverage

    def calc_median_coverage(self) -> float:
        # Method for calculating median coverage.

        median_cov: float
        if len(self._coverages) == 0:
            median_cov = None # no coverage values are available
        else:
            median_cov = median(self._coverages)
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

    # Number of termini of a contig
    num_contig_termini: int = 2
    # Total number of dead ends taking account of multiplicity
    total_dead_ends: int = 0

    i: int
    contig: Contig
    # Iterate over contigs
    for i, contig in enumerate(contig_collection):
        # Count overlaps associated with start
        start_is_not_dead: int = is_not_empty(tuple(
            filter(is_start_match, overlap_collection[i])
        ))
        # Count overlaps associated with end
        end_is_not_dead:   int = is_not_empty(tuple(
            filter(is_end_match, overlap_collection[i])
        ))

        # Calculate number of dead ends of the current contig
        num_dead_ends: int = num_contig_termini \
                             - start_is_not_dead - end_is_not_dead
        # Add to `total_dead_ends`
        total_dead_ends += num_dead_ends
    # end for

    # Total number of termini taking account of multiplicity
    total_termini: int = int(num_contig_termini * len(contig_collection))
    # Calculate the LQ coefficient
    lq_coef: float = (1 - total_dead_ends / total_termini) * 100.0

    return round(lq_coef, 2)
# end def calc_lq_coef


def calc_exp_genome_size(contig_collection: ContigCollection,
                         overlap_collection: OverlapCollection) -> int:
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;

    # In this variable, total length of overlapping regions will be stored
    total_overlap_len: int = 0

    # Iterate over contigs
    i: ContigIndex
    contig: Contig
    for i, contig in enumerate(contig_collection):

        # Get start-associated overlaps
        start_ovls: Sequence[Overlap] = tuple(
            filter(
                is_start_match,
                overlap_collection[i]
            )
        )
        # Get end-associated overlaps
        end_ovls: Sequence[Overlap] = tuple(
            filter(
                is_end_match,
                overlap_collection[i]
            )
        )

        # Function for `filter`.
        # Purpose: we won't consider contigs with index > i
        #   in order not to count an overlap twice.
        not_already_counted: Callable = lambda x: x.contig_j >= i

        ovls: Sequence[Overlap]
        for ovls in (start_ovls, end_ovls):

            ovls_to_add: Sequence[Overlap]

            if len(ovls) <= contig.multplty:
                # No extra overlaps. # We will just add lengths of overlaps to `total_overlap_len`.
                ovls_to_add = filter(not_already_counted, ovls)
            else:
                # Some extra overlaps discovered.
                # We will consider only M longest overlaps,
                #   where M is contig's multiplicity.
                ovls_to_add = filter(
                    not_already_counted,
                    sorted(                   # sort overlaps in a descending order
                        ovls, key=lambda x: -x.ovl_len
                    )[: int(contig.multplty)] # select M first (M the longest) overlaps
                )
            # end if

            # Add lengths of overlaps to `total_overlap_len`
            ovl: Overlap
            for ovl in ovls_to_add:
                total_overlap_len += ovl.ovl_len
            # end for
        # end for
    # end for

    # Calculate length of the genome, taking account of multiplicity of contigs.
    expected_genome_size: int = sum(                      # Sum products of contigs'
        map(                                              #   lengths and multiplicities:
            lambda x: int(x.length * round(x.multplty)),  # len1*multplty1 + len2*multplty2 + ...
            contig_collection
        )
    )\
    - total_overlap_len                            # Subtract total length of overlapping regions

    return expected_genome_size
# end def calc_exp_genome_size


def is_start_match(ovl: Overlap) -> bool:
    # Function returns True if overlap `ovl` is associated with start.
    return (ovl.terminus_i == START and ovl.terminus_j == END)\
           or (ovl.terminus_i == START and ovl.terminus_j == RCSTART)
           # or {ovl.terminus_i, ovl.terminus_j} == _START2START_MATCH
# end def def is_start_match


def is_end_match(ovl: Overlap) -> bool:
    # Function returns True if overlap `ovl` is associated with end.
    return (ovl.terminus_i == END and ovl.terminus_j == START)\
           or (ovl.terminus_i == END and ovl.terminus_j == RCEND)
           # or {ovl.terminus_i, ovl.terminus_j} == _END2END_MATCH
# end def def is_end_match


def is_not_empty(collection: Sequence) -> int:
    return int(len(collection) != 0)
# end def not_empty
