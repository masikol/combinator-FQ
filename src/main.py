# -*- coding: utf-8 -*-

from typing import Sequence, Dict, Any

import src.output as out
import src.contigs as cnt
import src.overlaps as ovl
from src.parse_args import parse_args
from src.filesystem import conf_prefix, make_outdir

def main(version, last_update_date):
    contigs_fpaths: Sequence[str]
    params: Dict[str, Any]

    contigs_fpaths, params = parse_args(version, last_update_date)

    fpath: str
    for fpath in contigs_fpaths:
        print('Processing file `{}`'.format(fpath))

        contig_collection: cnt.ContigCollection = cnt.get_contig_collection(fpath, params['a'])
        cnt.assign_multiplty(contig_collection)

        overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
            contig_collection, params['i'], params['a']
        )

        prefix: str = conf_prefix(fpath, params['o'])

        make_outdir(params['o'])
        out.write_adjacency_table(contig_collection, overlap_collection, params['o'], prefix)
        out.write_full_log(contig_collection, overlap_collection, params['o'], prefix)
        out.write_summary(contig_collection, overlap_collection, fpath, params['o'], prefix)
        print('-'*20)
    # end for
# end def main
