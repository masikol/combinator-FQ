# -*- coding: utf-8 -*-

from typing import Sequence, Dict, Any

import src.output as out
import src.contigs as cnt
import src.overlaps as ovl
from src.parse_args import parse_args
from src.filesystem import conf_prefix, make_outdir


def main(version: str, last_update_date: str) -> None:

    contigs_fpaths: Sequence[str] # paths to input files
    params: Dict[str, Any] # parameters of the program

    # Parse command line arguments
    contigs_fpaths, params = parse_args(version, last_update_date)

    # Report parameters of current run
    _report_parameters(params, version, last_update_date)

    # Iterate over input files and process them
    fpath: str
    for fpath in contigs_fpaths:
        print('Processing file `{}`'.format(fpath))

        # Create output dir
        make_outdir(params['o'])

        # Read contigs
        contig_collection: cnt.ContigCollection = cnt.get_contig_collection(fpath, params['a'])
        # Assign multiplicity to contigs
        cnt.assign_multiplty(contig_collection)

        # Detect adjacent contigs
        overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
            contig_collection, params['i'], params['a']
        )

        # Make prefix for current input file
        prefix: str = conf_prefix(fpath, params['o'])

        # Write output files
        # Write adjacency table
        out.write_adjacency_table(contig_collection, overlap_collection, params['o'], prefix)
        # Write full log
        out.write_full_log(contig_collection, overlap_collection, params['o'], prefix)
        # Write summary
        out.write_summary(contig_collection, overlap_collection, fpath, params['o'], prefix)

        print('-'*20)
    # end for
# end def main


def _report_parameters(params, version, last_update_date):
    print('{}. Version {}. {} edition.'.format('combinator-FQ', version, last_update_date))
    print('Parameters:')
    print(' - Minimum k: {} bp.'.format(params['i']))
    print(' - Maximum k: {} bp.'.format(params['a']))
    print(' - Output directory: `{}`.'.format(params['o']))
    print('-' * 20)
# end def _report_parameters
