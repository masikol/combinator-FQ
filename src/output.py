# -*- coding: utf-8 -*-

import os
import sys
from typing import TextIO, Callable

import src.statistics as sts
from src.contigs import ContigCollection
from src.overlaps import OverlapCollection

def write_summary(contig_collection: ContigCollection,
                  overlap_collection: OverlapCollection,
                  infpath: str, outdpath: str, out_prefix: str) -> None:

    summary_fpath: str = os.path.join(
        outdpath,
        "{}{}combinator_summary_FQ.txt"\
            .format(out_prefix, '' if out_prefix == '' else '_')
    )

    out_str: str
    outfile: TextIO
    with open(summary_fpath, 'w') as outfile:

        outfile.write("File: '{}'\n\n".format(infpath))

        # Summary with some statistics:
        _double_write(' === Summary ===\n', outfile)

        # Number of contigs processed:
        out_str = '{} contigs were processed.'.format(len(contig_collection))
        _double_write(out_str, outfile)

        # Sum of contigs' lengths
        out_str = 'Sum of contig lengths: {} b.p.'\
            .format(sts.calc_sum_contig_lengths(contig_collection))
        _double_write(out_str, outfile)

        # Expected length of the genome
        out_str = 'Expected length of the genome: {} b.p.'.format(-1) # STUB !!
        _double_write(out_str, outfile)

        # Calculate coverage statistics
        cov_calc = sts.CoverageCalculator(contig_collection)

        # Min coverage
        min_coverage: float = cov_calc.get_min_coverage()
        out_str = 'Min coverage: {}'\
            .format(min_coverage if not min_coverage is None else 'NA')
        _double_write(out_str, outfile)

        # Max coverage
        max_coverage: float = cov_calc.get_max_coverage()
        out_str = 'Max coverage: {}'\
            .format(max_coverage if not max_coverage is None else 'NA')
        _double_write(out_str, outfile)

        # Mean coverage
        mean_coverage: float = cov_calc.calc_mean_coverage()
        out_str = 'Mean coverage: {}'\
            .format(mean_coverage if not mean_coverage is None else 'NA')
        _double_write(out_str, outfile)

        # LQ coefficient
        out_str = 'LQ-coefficient: {}'\
            .format(sts.calc_lq_coef(contig_collection, overlap_collection))
        _double_write(out_str, outfile)
    # end with
# end def write_summary


def _double_write(outstr: str, outfile: TextIO) -> None:
    print_func: Callable[[str], int]
    for print_func in (sys.stdout.write, outfile.write):
        print_func('{}\n'.format(outstr))
    # end for
    sys.stdout.flush()
# end def double_write