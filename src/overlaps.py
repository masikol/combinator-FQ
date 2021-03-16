# -*- coding: utf-8 -*-

from src.find_overlap import find_overlap_s2s, find_overlap_e2s, find_overlap_e2e


START, RCSTART, END, RCEND = range(4)

_KEY2LETTER_MAP = {
    0: 'S',
    1: 'rc_S',
    2: 'E',
    3: 'rc_E'
}

_KEY2WORD_MAP = {
    0: 'start',
    1: 'rc-start',
    2: 'end',
    3: 'rc-end'
}


class Overlap:

    def __init__(self, contig1, terminus1, contig2, terminus2, ovl_len):
        self.contig1 = contig1
        self.terminus1 = terminus1
        self.contig2 = contig2
        self.terminus2 = terminus2
        self.ovl_len = ovl_len
    # end def __init__

    def __repr__(self):
        return '<{}-{}; {}-{}; len={}>'\
        .format(self.contig1, self.terminus1,
        self.contig2, self.terminus2, self.ovl_len)
    # end def __repr__
# end class Overlap


class OverlapCollection:

    def __init__(self):
        self._collection = dict()
    # end def

    def __getitem__(self, position):
        try:
            item_to_return = self._collection[position]
        except KeyError:
            item_to_return = tuple()
        finally:
            return item_to_return
        # end try
    # end def __getitem__

    def add_overlap(self, key, overlap):
        try:
            self._collection[key].append(overlap)
        except KeyError:
            self._collection[key] = [overlap]
        # end try
    # end def add_overlap
# end class OverlapCollection


def get_overlaps_str_for_table(overlap_collection, contig_collection, key):

    if len(overlap_collection[key]) == 0:
        return '-'
    else:
        overlaps = overlap_collection[key]
        match_strings = list()

        for ovl in overlaps:
            if ovl.contig1 != ovl.contig2:
                letter1 = _KEY2LETTER_MAP[ovl.terminus1]
                letter2 = _KEY2LETTER_MAP[ovl.terminus2]
                match_strings.append('[{}={}({}); ovl={}]'\
                    .format(letter1, letter2, contig_collection[ovl.contig2].name, ovl.ovl_len))
            else:
                match_strings.append('[Circle; ovl={}]'.format(ovl.ovl_len))
            # end if
        # end for

        return ' '.join(match_strings)
    # end if
# end def get_overlaps_str_for_table


def get_overlaps_str_for_log(overlap_collection, contig_collection, key):
    if len(overlap_collection[key]) == 0:
        return ''
    else:
        overlaps = overlap_collection[key]
        match_strings = list()

        for ovl in overlaps:
            if ovl.contig1 != ovl.contig2:
                word1 = _KEY2WORD_MAP[ovl.terminus1]
                word2 = _KEY2WORD_MAP[ovl.terminus2]
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


def detect_adjacent_contigs(contig_collection, mink, maxk):

    num_contigs = len(contig_collection)
    overlap_collection = OverlapCollection()

    for i in range(num_contigs):

        # Omit contigs shorter that 'mink'
        if contig_collection[i].length <= mink:
            print('\r{}/{}'.format(i+1, num_contigs), end='')
            continue
        # end if

        # === Compare start of the current contig to end of the current contig ===
        ovl_len = find_overlap_e2s(contig_collection[i].end,
                                   contig_collection[i].start,
                                   mink, maxk)
        if ovl_len != 0 and contig_collection[i].length != ovl_len:
            overlap_collection.add_overlap(i, Overlap(i, END, i, START, ovl_len))
            overlap_collection.add_overlap(i, Overlap(i, START, i, END, ovl_len))
        # end if

        # === Compare start of the current conitg to rc-end of the current contig ===
        ovl_len = find_overlap_s2s(contig_collection[i].start,
                                   contig_collection[i].rcend,
                                   mink, maxk)
        if ovl_len != 0:
            overlap_collection.add_overlap(i, Overlap(i, START, i, RCEND, ovl_len))
            overlap_collection.add_overlap(i, Overlap(i, RCEND, i, START, ovl_len))
        # end if


        # |=== Compare i-th contig to contigs from i+1 to N ===|
        #   in order not to compare pairs of contigs more than one time
        for j in range(i+1, num_contigs):

            # === Compare i-th start to j-th end ===
            ovl_len = find_overlap_e2s(contig_collection[j].end,
                                       contig_collection[i].start,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, END, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, END, i, START, ovl_len))
            # end if

            # === Compare i-th end to j-th start ===
            ovl_len = find_overlap_e2s(contig_collection[i].end,
                                       contig_collection[j].start,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, START, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, START, i, END, ovl_len))
            # end if

            # === Compare i-th start to reverse-complement j-th start ===
            ovl_len = find_overlap_e2s(contig_collection[j].rcstart,
                                       contig_collection[i].start,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCSTART, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCSTART, i, START, ovl_len))
            # end if

            # === Compare i-th end to reverse-complement j-th end ===
            ovl_len = find_overlap_e2s(contig_collection[i].end,
                                       contig_collection[j].rcend,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCEND, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCEND, i, END, ovl_len))
            # end if

            # === Compare i-th start to j-th start ===
            ovl_len = find_overlap_s2s(contig_collection[i].start,
                                       contig_collection[j].start,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, START, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, START, i, START, ovl_len))
            # end if

            # === Compare i-th end to j-th end ===
            ovl_len = find_overlap_e2e(contig_collection[i].end,
                                       contig_collection[j].end,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, END, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, END, i, END, ovl_len))
            # end if

            # === Compare i-th start to reverse-complement j-th end ===
            ovl_len = find_overlap_s2s(contig_collection[i].start,
                                       contig_collection[j].rcend,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCEND, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCEND, i, START, ovl_len))
            # end if

            # === Compare i-th end to reverse-complement j-th start ===
            ovl_len = find_overlap_e2e(contig_collection[i].end,
                                       contig_collection[j].rcstart,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCSTART, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCSTART, i, END, ovl_len))
            # end if
        # end for

        print('\r{}/{}'.format(i+1, num_contigs), end='')
    # end for
    print()

    return overlap_collection
# end def detect_adjacent_contigs
