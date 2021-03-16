# -*- coding: utf-8 -*-

from typing import NewType, Dict, List
from src.contigs import Contig, ContigCollection, ContigIndex
from src.find_overlap import find_overlap_s2s, find_overlap_e2s, find_overlap_e2e


Terminus = NewType('Terminus', int)


START:   Terminus = 0
RCSTART: Terminus = 1
END:     Terminus = 2
RCEND:   Terminus = 3

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


class Overlap:

    def __init__(self,
                 contig1: ContigIndex, terminus1: Terminus,
                 contig2: ContigIndex, terminus2: Terminus,
                 ovl_len: int):
        self.contig1 = contig1
        self.terminus1 = terminus1
        self.contig2 = contig2
        self.terminus2 = terminus2
        self.ovl_len: int = ovl_len
    # end def __init__

    # def __repr__(self):
    #     return '<{}-{}; {}-{}; len={}>'\
    #     .format(self.contig1, self.terminus1,
    #     self.contig2, self.terminus2, self.ovl_len)
    # # end def __repr__
# end class Overlap


class OverlapCollection:

    def __init__(self) -> None:
        self._collection: Dict[List[Overlap]] = dict()
    # end def

    def __getitem__(self, key: ContigIndex) -> List[Overlap]:
        try:
            item_to_return = self._collection[key]
        except KeyError:
            # Return tuple just for efficiecy,
            #   though, strictly speaking, it should have been a list.
            item_to_return = tuple()
        # end try
        return item_to_return
    # end def __getitem__

    def add_overlap(self, key: ContigIndex, overlap: Overlap) -> None:
        try:
            self._collection[key].append(overlap)
        except KeyError:
            self._collection[key] = [overlap]
        # end try
    # end def add_overlap
# end class OverlapCollection


def get_overlaps_str_for_table(overlap_collection: OverlapCollection,
                               contig_collection:  ContigCollection,
                               key: ContigIndex) -> str:

    if len(overlap_collection[key]) == 0:
        return '-'
    else:
        overlaps: List[Overlap] = overlap_collection[key]
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
# end def get_overlaps_str_for_table


def get_overlaps_str_for_log(overlap_collection: OverlapCollection,
                             contig_collection:  ContigCollection,
                             key: ContigIndex) -> str:
    if len(overlap_collection[key]) == 0:
        return ''
    else:
        overlaps: List[Overlap] = overlap_collection[key]
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


def detect_adjacent_contigs(contig_collection: ContigCollection,
                            mink: int, maxk: int) -> OverlapCollection:

    num_contigs: int = len(contig_collection)
    overlap_collection: OverlapCollection = OverlapCollection()

    i: ContigIndex
    for i in range(num_contigs):

        # Omit contigs shorter that 'mink'
        if contig_collection[i].length <= mink:
            print('\r{}/{}'.format(i+1, num_contigs), end='')
            continue
        # end if

        # === Compare start of the current contig to end of the current contig ===
        ovl_len: int = find_overlap_e2s(contig_collection[i].end,
                                        contig_collection[i].start,
                                        mink, maxk)
        if not ovl_len in (0, contig_collection[i].length):
            overlap_collection.add_overlap(i, Overlap(i, END, i, START, ovl_len))
            overlap_collection.add_overlap(i, Overlap(i, START, i, END, ovl_len))
        # end if

        # === Compare start of the current conitg to rc-end of the current contig ===
        ovl_len: int = find_overlap_s2s(contig_collection[i].start,
                                        contig_collection[i].rcend,
                                        mink, maxk)
        if ovl_len != 0:
            overlap_collection.add_overlap(i, Overlap(i, START, i, RCEND, ovl_len))
            overlap_collection.add_overlap(i, Overlap(i, RCEND, i, START, ovl_len))
        # end if


        # |=== Compare i-th contig to contigs from i+1 to N ===|
        #   in order not to compare pairs of contigs more than one time
        j: ContigIndex
        for j in range(i+1, num_contigs):

            # === Compare i-th start to j-th end ===
            ovl_len: int = find_overlap_e2s(contig_collection[j].end,
                                            contig_collection[i].start,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, END, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, END, i, START, ovl_len))
            # end if

            # === Compare i-th end to j-th start ===
            ovl_len: int = find_overlap_e2s(contig_collection[i].end,
                                            contig_collection[j].start,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, START, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, START, i, END, ovl_len))
            # end if

            # === Compare i-th start to reverse-complement j-th start ===
            ovl_len: int = find_overlap_e2s(contig_collection[j].rcstart,
                                            contig_collection[i].start,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCSTART, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCSTART, i, START, ovl_len))
            # end if

            # === Compare i-th end to reverse-complement j-th end ===
            ovl_len: int = find_overlap_e2s(contig_collection[i].end,
                                            contig_collection[j].rcend,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCEND, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCEND, i, END, ovl_len))
            # end if

            # === Compare i-th start to j-th start ===
            ovl_len: int = find_overlap_s2s(contig_collection[i].start,
                                            contig_collection[j].start,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, START, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, START, i, START, ovl_len))
            # end if

            # === Compare i-th end to j-th end ===
            ovl_len: int = find_overlap_e2e(contig_collection[i].end,
                                            contig_collection[j].end,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, END, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, END, i, END, ovl_len))
            # end if

            # === Compare i-th start to reverse-complement j-th end ===
            ovl_len: int = find_overlap_s2s(contig_collection[i].start,
                                            contig_collection[j].rcend,
                                            mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCEND, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, RCEND, i, START, ovl_len))
            # end if

            # === Compare i-th end to reverse-complement j-th start ===
            ovl_len: int = find_overlap_e2e(contig_collection[i].end,
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
