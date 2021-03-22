# -*- coding: utf-8 -*-

import os
import sys
from typing import TextIO, Callable, Dict, Collection, List, MutableSequence

from src.platform import platf_depend_exit
import src.combinator_statistics as sts
from src.contigs import Contig, ContigCollection, ContigIndex
from src.overlaps import Overlap, OverlapCollection
from src.overlaps import Terminus, START, RCSTART, END, RCEND


# Dictionary maps `Terminus` to it's "letter" representation
#   for adjacency table.
_KEY2LETTER_MAP: Dict[Terminus, str] = {
    START: 'S',
    RCSTART: 'rc_S',
    END: 'E',
    RCEND: 'rc_E'
}

# Dictionary maps `Terminus` to it's "word" representation
#   for full log.
_KEY2WORD_MAP: Dict[Terminus, str] = {
    START: 'start',
    RCSTART: 'rc-start',
    END: 'end',
    RCEND: 'rc-end'
}


def write_summary(contig_collection: ContigCollection,
                  overlap_collection: OverlapCollection,
                  infpath: str, outdpath: str, out_prefix: str) -> None:
    # Function writes summary to summary file.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;
    # :param infpath: path to input file (it will be mentioned in summary);
    # :param outdpath: path to output directory;
    # :param out_prefix: prefix for current output files;

    # Make path to summary file
    summary_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_summary_FQ.txt'\
            .format(out_prefix)
    )

    print('Writing summary to `{}`'.format(summary_fpath))

    # Proceed
    wrk_str: str
    outfile: TextIO
    with open(summary_fpath, 'w') as outfile:

        # Path to input file
        outfile.write('Input file: `{}`\n\n'.format(infpath))

        # Summary with some statistics:
        _double_write(' === Summary ===', outfile)

        # Number of contigs processed:
        wrk_str = '{} contigs were processed.'.format(len(contig_collection))
        _double_write(wrk_str, outfile)

        # Sum of contigs' lengths
        wrk_str = 'Sum of contig lengths: {} bp'\
            .format(sts.calc_sum_contig_lengths(contig_collection))
        _double_write(wrk_str, outfile)

        # Expected length of the genome
        wrk_str = 'Expected length of the genome: {} bp'.format(-1) # STUB !!
        _double_write(wrk_str, outfile)

        # Calculate coverage statistics
        cov_calc = sts.CoverageCalculator(contig_collection)

        # Min coverage
        min_coverage: float = cov_calc.get_min_coverage()
        wrk_str = 'Min coverage: {}'\
            .format(min_coverage if not min_coverage is None else 'NA')
        _double_write(wrk_str, outfile)

        # Max coverage
        max_coverage: float = cov_calc.get_max_coverage()
        wrk_str = 'Max coverage: {}'\
            .format(max_coverage if not max_coverage is None else 'NA')
        _double_write(wrk_str, outfile)

        # Mean coverage
        mean_coverage: float = cov_calc.calc_mean_coverage()
        wrk_str = 'Mean coverage: {}'\
            .format(mean_coverage if not mean_coverage is None else 'NA')
        _double_write(wrk_str, outfile)

        # Median coverage
        median_coverage: float = cov_calc.calc_median_coverage()
        wrk_str = 'Median coverage: {}'\
            .format(median_coverage if not median_coverage is None else 'NA')
        _double_write(wrk_str, outfile)

        # LQ coefficient
        wrk_str = 'LQ-coefficient: {}'\
            .format(sts.calc_lq_coef(contig_collection, overlap_collection))
        _double_write(wrk_str, outfile)
    # end with
# end def write_summary


def write_adjacency_table(contig_collection: ContigCollection,
                          overlap_collection: OverlapCollection,
                          outdpath: str, out_prefix: str) -> None:
    # Function writes adjacency table to output TSV file.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;
    # :param outdpath: path to output directory;
    # :param out_prefix: prefix for current output files;

    # Make path to output TSV file
    adj_table_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'\
            .format(out_prefix)
    )

    print('Writing adjacency table to `{}`'.format(adj_table_fpath))

    # Proceed
    outfile: TextIO
    with open(adj_table_fpath, 'w') as outfile:

        # Write head of the table:
        outfile.write('\t'.join([
            '#',
            'Contig name',
            'Length',
            'Coverage',
            'GC(%)',
            'Multiplicity',
            'Annotation',
            'Start',
            'End',
        ]))
        outfile.write('\n')

        wrk_str: str

        # Iterate over contigs and write their properties
        i: ContigIndex
        contig: Contig
        for i, contig in enumerate(contig_collection):

            # Ordinal number and name
            outfile.write('{}\t'.format(i + 1))
            outfile.write('{}\t'.format(contig.name))

            # Length
            outfile.write('{}\t'.format(contig.length))

            # Coverage
            wrk_str = '-' if contig.cov is None else str(contig.cov)
            outfile.write('{}\t'.format(wrk_str))

            # GC content of the contig
            outfile.write(str(contig.gc_content) + '\t')

            # Multiplicity
            wrk_str = '-' if contig.cov is None else str(contig.multplty)
            outfile.write('{}\t'.format(wrk_str))

            # Empty column for annotation
            outfile.write('\t')

            # Information about discovered adjacency
            # "Start" column
            wrk_str = _get_overlaps_str_for_table(overlap_collection,
                                                  contig_collection,
                                                  i, 's')
            outfile.write('{}\t'.format(wrk_str))

            # "End" column
            wrk_str = _get_overlaps_str_for_table(overlap_collection,
                                                  contig_collection,
                                                  i, 'e')
            outfile.write('{}'.format(wrk_str))

            outfile.write('\n')
        # end for
    # end with
# end def write_adjacency_table


def write_full_log(contig_collection: ContigCollection,
                   overlap_collection: OverlapCollection,
                   outdpath: str, out_prefix: str) -> None:
    # Function writes full matching log (not only adjacency-associated matches)
    #   to "full-log" file.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;
    # :param outdpath: path to output directory;
    # :param out_prefix: prefix for current output files;

    # Make path to full log file
    log_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_full_matching_log.txt'\
            .format(out_prefix)
    )

    print('Writing full matching log to `{}`'.format(log_fpath))

    # Proceed
    outfile: TextIO
    with open(log_fpath, 'w') as outfile:

        wrk_str: str

        # Write information about discovered adjacency
        i: ContigIndex
        for i, _ in enumerate(contig_collection):
            # Write what matches start of current contig
            wrk_str = _get_overlaps_str_for_log(overlap_collection,
                                                contig_collection, i)
            if not wrk_str is None:
                outfile.write('{}\n'.format(wrk_str))
            # end if
        # end for
    # end with
# end def write_full_log


def _double_write(outstr: str, outfile: TextIO) -> None:
    # Function for writing and printing.
    # "Double" means that it writes identical data both to
    #   `outfile` and to stdout.
    #
    # :param outstr: string to be printed/written;
    # :param outfile: file-like instance of output file to write in;
    print_func: Callable[[str], int]
    for print_func in (sys.stdout.write, outfile.write):
        print_func('{}\n'.format(outstr))
    # end for
    sys.stdout.flush()
# end def double_write


def _get_start_matches(overlaps: MutableSequence[Overlap]) -> Collection[Overlap]:
    # Function selects "start-associated" overlaps from a collection of overlaps.
    return tuple(
        filter(sts.is_start_match, overlaps)
    )
# end def _get_start_matches


def _get_end_matches(overlaps: MutableSequence[Overlap]) -> Collection[Overlap]:
    # Function selects "end-associated" overlaps from a collection of overlaps.
    return tuple(
        filter(sts.is_end_match, overlaps)
    )
# end def _get_end_matches


def _select_get_matches(term: str) -> Callable[[MutableSequence[Overlap]], Collection[Overlap]]:
    # Function returns function depending on `term` (terminus) parameter.
    # If `term` is 's', it returns `_get_start_matches` function.
    # If `term` is 'e', it returns `_get_end_matches` function.
    #
    # :param term: string 's' (start) or 'e' (end);

    if term == 's':
        get_matches = _get_start_matches
    elif term == 'e':
        get_matches = _get_end_matches
    else:
        print('Fatal error: invalid value passed to function \
`_get_overlaps_str_for_table` with argument `term`: `{}`'.format(term))
        print('Please, contact the developer.')
        platf_depend_exit(1)
    # end if

    return get_matches
# end def _select_get_matches


def _get_overlaps_str_for_table(overlap_collection: OverlapCollection,
                                contig_collection:  ContigCollection,
                                key: ContigIndex, term: str) -> str:
    # Function extracts overlaps of `key` contigs associated with `term` terminus,
    #   converts this collection of `src.overlaps.Overlap` to string representation
    #   for adjacency table.
    # After this, the fucntion returns this string.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;
    # :param key: key (index) of contig;
    # :param term: "terminus" -- 's' or 'e';

    # Select function for obtaining `term`-associated overlaps.
    get_matches: Callable[[List[Overlap]], Collection[Overlap]]
    get_matches = _select_get_matches(term)

    # Extract overlaps associated with `term` terminus for `key` contig
    overlaps: Collection[Overlap] = get_matches(overlap_collection[key])

    # Convert `Overlap` instances to string representation
    if len(overlaps) == 0:
        return '-' # no proper overlaps found
    else:
        match_strings: List[str] = list() # list for formatted strings

        ovl: Overlap
        for ovl in overlaps:
            # If contig does not match itself
            if ovl.contig_i != ovl.contig_j:
                # Letter for the first contig of the overlap
                letter1: str = _KEY2LETTER_MAP[ovl.terminus_i]
                # Letter for the second contig of the overlap
                letter2: str = _KEY2LETTER_MAP[ovl.terminus_j]
                # Convert and append
                match_strings.append('[{}={}({}); ovl={}]'\
                    .format(letter1, letter2, contig_collection[ovl.contig_j].name,
                            ovl.ovl_len))
            # If contig matches itself (it is circular)
            else:
                match_strings.append('[Circle; ovl={}]'.format(ovl.ovl_len))
            # end if
        # end for

        # Separate formatted strings with spaces
        return ' '.join(match_strings)
    # end if
# end def _get_overlaps_str_for_table


def _get_overlaps_str_for_log(overlap_collection: OverlapCollection,
                              contig_collection:  ContigCollection,
                              key: ContigIndex) -> str:
    # Function extracts overlaps of `key` contigs associated with `term` terminus,
    #   converts this collection of `src.overlaps.Overlap` to string representation
    #   for full log.
    # After this, the fucntion returns this string.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param overlap_collection: instance of OverlapCollection returned by
    #   `src.overlaps.detect_adjacent_contigs` function;
    # :param key: key (index) of contig;

    # Extract overlaps for current contig
    overlaps: List[Overlap] = overlap_collection[key]

    if len(overlaps) == 0:
        return None # no proper overlaps found
    else:
        match_strings: List[str] = list() # list for formatted strings

        ovl: Overlap
        for ovl in overlaps:
            # If contig does not match itself
            if ovl.contig_i != ovl.contig_j:

                # Word for the first contig of the overlap
                word1: str = _KEY2WORD_MAP[ovl.terminus_i]
                # Word for the second contig of the overlap
                word2: str = _KEY2WORD_MAP[ovl.terminus_j]

                # Convert and append
                match_strings.append('{}: {} matches {} of {} with overlap of {} bp'\
                    .format(contig_collection[key].name, word1, word2,
                        contig_collection[ovl.contig_j].name, ovl.ovl_len))
            else:
                # Contig is circular
                if ovl.terminus_i == END and ovl.terminus_j == START:
                    match_strings.append('{}: contig is circular with overlap of {} bp'\
                        .format(contig_collection[key].name, ovl.ovl_len))
                # Start of contig matches it's own reverse-complement end
                elif ovl.terminus_i == START and ovl.terminus_j == RCEND:
                    match_strings.append('{}: start is identical to it\'s own rc-end with overlap of {} bp'\
                        .format(contig_collection[key].name, ovl.ovl_len))
                # end if
            # end if
        # end for

        # Separate formatted strings with new line chars
        return '\n'.join(match_strings)
    # end if
# end def
