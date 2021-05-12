# -*- coding: utf-8 -*-

from typing import MutableSequence

import src.combinator_statistics as sts
from src.contigs import ContigCollection, ContigIndex
from src.overlaps import OverlapCollection, Overlap


def assign_multiplty(contig_collection: ContigCollection, overlap_collection: OverlapCollection) -> None:
    # Function assigns multiplicity (copies of this contig in the genome) to contigs.
    # :param contig_collection: instance of `ContigCollection`
    #   returned by function `get_contig_collection`;
    # Function modifies `contig_collection` argument passed to it.

    # Coverage of 1-st contig can be zero.
    # In this case we cannot calculate multiplicity of contigs based on coverage.
    first_cov_is_valid: bool = True

    # Set `first_cov_is_valid` to False if the first contig has no coverage.
    if contig_collection[0].cov is None:
        first_cov_is_valid = False
    # end if

    # Set `first_cov_is_valid` to False if the first contig has zero coverage.
    # And report it.
    if first_cov_is_valid and contig_collection[0].cov < 1e-6:
        # Coverage of 1-st contig is zero
        print('\n`{}` has zero coverage (less than 1e-6 actually).'\
            .format(contig_collection[0].name))
        print('Multiplicity of contigs will be calculated based on overlaps instead of coverage.\n')
        first_cov_is_valid = False
    # end if

    # Calculate multiplicity of contigs:
    i: ContigIndex
    for i in range(len(contig_collection)):

        # Validate coverage of the first contig
        calc_multiplty_by_div: bool = first_cov_is_valid

        # Validate coverage of the current contig
        calc_multiplty_by_div = calc_multiplty_by_div and not contig_collection[i].cov is None

        if calc_multiplty_by_div:
            # Calculate multiplicity based on coverage
            contig_collection[i].multplty = _calc_multiplty_by_coverage(
                contig_collection[i].cov,
                contig_collection[0].cov
            )
        else:
            # Calculate multiplicity based on overlaps
            contig_collection[i].multplty = _calc_multiplty_by_overlaps(
                overlap_collection[i]
            )
        # end if
    # end for
# end def assign_multiplty


def _calc_multiplty_by_coverage(curr_coverage: float, first_contig_coverage: float) -> float:
    # Function for calculating multiplicity of a given contig
    #   based on it's coverage (and also on coverage of the first contig).
    # :param curr_coverage: coverage of current contig;
    # :param curr_coverage: coverage of the first contig;

    multiplicity: float = curr_coverage / first_contig_coverage
    return max(multiplicity, 1.0)
# end def _calc_multiplty_by_coverage


def _calc_multiplty_by_overlaps(ovl_list: MutableSequence[Overlap]) -> float:
    # Function for calculating multiplicity of a given contig
    #   based on number of overlaps of this contig.
    # :param ovl_list: list of overlaps of current contig;

    # Count overlaps associated with start
    num_start_matches: int = len(tuple(
        filter(sts.is_start_match, ovl_list)
    ))
    # Count overlaps associated with end
    num_end_matches:   int = len(tuple(
        filter(sts.is_end_match, ovl_list)
    ))

    # Obtain multiplicity based on number of overlaps
    multiplicity: float = max(
        1.0,
        min(num_start_matches, num_end_matches)
    )

    return multiplicity
# end def _calc_multiplty_by_overlaps
