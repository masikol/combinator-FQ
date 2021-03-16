# -*- coding: utf-8 -*-

import os
import sys
from typing import TextIO, Callable, Dict, Collection, List

import src.statistics as sts
from src.contigs import Contig, ContigCollection, ContigIndex
from src.overlaps import Overlap, OverlapCollection
from src.overlaps import Terminus, START, RCSTART, END, RCEND


_KEY2LETTER_MAP: Dict[Terminus, str] = {
    START: 'S',
    RCSTART: 'rc_S',
    END: 'E',
    RCEND: 'rc_E'
}

_KEY2WORD_MAP: Dict[Terminus, str] = {
    START: 'start',
    RCSTART: 'rc-start',
    END: 'end',
    RCEND: 'rc-end'
}


def write_summary(contig_collection: ContigCollection,
                  overlap_collection: OverlapCollection,
                  infpath: str, outdpath: str, out_prefix: str) -> None:

    summary_fpath: str = os.path.join(
        outdpath,
        '{}{}combinator_summary_FQ.txt'\
            .format(out_prefix, '' if out_prefix == '' else '_')
    )

    print('Writing summary to `{}`'.format(summary_fpath))

    out_str: str
    outfile: TextIO
    with open(summary_fpath, 'w') as outfile:

        outfile.write('Input file: `{}`\n\n'.format(infpath))

        # Summary with some statistics:
        _double_write(' === Summary ===', outfile)

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


def write_adjacency_table(contig_collection: ContigCollection,
                          overlap_collection: OverlapCollection,
                          outdpath: str, out_prefix: str) -> None:

    adj_table_fpath: str = os.path.join(
        outdpath,
        '{}{}combinator_adjacent_contigs.tsv'\
            .format(out_prefix, '' if out_prefix == '' else '_')
    )

    print('Writing adjacency table to `{}`'.format(adj_table_fpath))

    outfile: TextIO
    with open(adj_table_fpath, 'w') as outfile:

        # Write head of the table:
        outfile.write("#\tContig name\tLength\tCoverage\tGC(%)\tMultiplicity\tAnnotation\tStart\tEnd\n")

        buff_out_str: str

        i: ContigIndex
        contig: Contig
        for i, contig in enumerate(contig_collection):

            # Write ordinal number and name
            outfile.write('{}\t'.format(i + 1))
            outfile.write('{}\t'.format(contig.name))

            # Write length
            outfile.write('{}\t'.format(contig.length))

            # Write coverage
            buff_out_str = '-' if contig.cov is None else str(contig.cov)
            outfile.write('{}\t'.format(buff_out_str))

            # Write GC content of the contig
            outfile.write(str(contig.gc_content) + '\t')

            # Write multiplicity
            buff_out_str = '-' if contig.cov is None else str(contig.multplty)
            outfile.write('{}\t'.format(buff_out_str))

            # Write empty column for annotation
            outfile.write('\t')

            # Write information about discovered adjacency
            # Write "Start" column
            buff_out_str = _get_overlaps_str_for_table(overlap_collection,
                                                       contig_collection,
                                                       i, 's')
            outfile.write('{}\t'.format(buff_out_str))

            # Write "End" column
            buff_out_str = _get_overlaps_str_for_table(overlap_collection,
                                                       contig_collection,
                                                       i, 'e')
            outfile.write('{}'.format(buff_out_str))

            outfile.write('\n')
        # end for
    # end with
# end def write_adjacency_table


def write_full_log(contig_collection: ContigCollection,
                   overlap_collection: OverlapCollection,
                   outdpath: str, out_prefix: str) -> None:

    log_fpath: str = os.path.join(
        outdpath,
        '{}{}combinator_full_matching_log.txt'\
            .format(out_prefix, '' if out_prefix == '' else '_')
    )

    print('Writing full matching log to `{}`'.format(log_fpath))

    outfile: TextIO
    with open(log_fpath, 'w') as outfile:

        buff_out_str: str

        # Write information about discovered adjacency
        i: ContigIndex
        contig: Contig
        for i, contig in enumerate(contig_collection):
            # Write what matches start of current contig
            buff_out_str = _get_overlaps_str_for_log(overlap_collection,
                                                     contig_collection, i)
            if buff_out_str != '':
                outfile.write('{}\n'.format(buff_out_str))
            # end if
        # end for
    # end with
# end def write_full_log


def _double_write(outstr: str, outfile: TextIO) -> None:
    print_func: Callable[[str], int]
    for print_func in (sys.stdout.write, outfile.write):
        print_func('{}\n'.format(outstr))
    # end for
    sys.stdout.flush()
# end def double_write


def _get_start_matches(overlaps: List[Overlap]) -> Collection[Overlap]:
    return tuple(
        filter(sts.is_start_match, overlaps)
    )
# end def _get_start_matches


def _get_end_matches(overlaps: List[Overlap]) -> Collection[Overlap]:
    return tuple(
        filter(sts.is_end_match, overlaps)
    )
# end def _get_end_matches


def _select_get_matches(term: str) -> Callable[[List[Overlap]], Collection[Overlap]]:
    if term == 's':
        get_matches = _get_start_matches
    elif term == 'e':
        get_matches = _get_end_matches
    else:
        raise ValueError
    # end if
    return get_matches
# end def _select_get_matches


def _get_overlaps_str_for_table(overlap_collection: OverlapCollection,
                                contig_collection:  ContigCollection,
                                key: ContigIndex, term: str) -> str:

    get_matches: Callable[[List[Overlap]], Collection[Overlap]]
    try:
        get_matches = _select_get_matches(term)
    except ValueError:
        raise ValueError('Invalid value passed to function \
`_get_overlaps_str_for_table` with argument `term`: `{}`'.format(term))
    # end try

    overlaps: Collection[Overlap] = get_matches(overlap_collection[key])

    if len(overlaps) == 0:
        return '-'
    else:
        match_strings: List[str] = list()

        ovl: Overlap
        for ovl in overlaps:
            if ovl.contig1 != ovl.contig2:
                letter1: str = _KEY2LETTER_MAP[ovl.terminus1]
                letter2: str = _KEY2LETTER_MAP[ovl.terminus2]
                match_strings.append('[{}={}({}); ovl={}]'\
                    .format(letter1, letter2, contig_collection[ovl.contig2].name, ovl.ovl_len))
            else:
                match_strings.append('[Circle; ovl={}]'.format(ovl.ovl_len))
            # end if
        # end for

        return ' '.join(match_strings)
    # end if
# end def _get_overlaps_str_for_table


def _get_overlaps_str_for_log(overlap_collection: OverlapCollection,
                              contig_collection:  ContigCollection,
                              key: ContigIndex) -> str:

    overlaps: List[Overlap] = overlap_collection[key]

    if len(overlaps) == 0:
        return ''
    else:
        match_strings: List[str] = list()

        ovl: Overlap
        for ovl in overlaps:

            if ovl.contig1 != ovl.contig2:

                word1: str = _KEY2WORD_MAP[ovl.terminus1]
                word2: str = _KEY2WORD_MAP[ovl.terminus2]

                match_strings.append('{}: {} matches {} of {} with overlap of {} b.p.'\
                    .format(contig_collection[key].name, word1, word2,
                        contig_collection[ovl.contig2].name, ovl.ovl_len))
            else:

                if ovl.terminus1 == END and ovl.terminus2 == START:
                    match_strings.append('{}: contig is circular with overlap of {} b.p.'\
                        .format(contig_collection[key].name, ovl.ovl_len))

                elif ovl.terminus1 == START and ovl.terminus2 == RCEND:
                    match_strings.append('{}: start is identical to it\'s own rc-end with overlap of {} b.p.'\
                        .format(contig_collection[key].name, ovl.ovl_len))
                # end if
            # end if
        # end for

        return '\n'.join(match_strings)
    # end if
# end def
