# -*- encoding: utf-8 -*-

from typing import NewType, Dict, List

from src.contigs import ContigCollection, ContigIndex
from src.find_overlap import find_overlap_s2s, find_overlap_e2s, find_overlap_e2e


# Following termini are defined in the program:
#   start, reverse-complement start, end, reverse-complement end.
Terminus = NewType('Terminus', int)


# Assign values for defined termini
START:   Terminus = 0
RCSTART: Terminus = 1
END:     Terminus = 2
RCEND:   Terminus = 3


class Overlap:
    # Class represents overlap between two contigs.

    def __init__(self,
                 contig_i: ContigIndex, terminus_i: Terminus,
                 contig_j: ContigIndex, terminus_j: Terminus,
                 ovl_len: int) -> None:
        # :param contig_i: index (key) of the 1-st contig;
        # :param contig_j: index (key) of the 2-nd contig;
        # :param terminus_i: terminus (of the 1-st contig) involved in the overlap;
        # :param terminus_j: terminus (of the 2-nd contig) involved in the overlap;
        # :param ovl_len: length of the overlap;
        self.contig_i = contig_i
        self.terminus_i = terminus_i
        self.contig_j = contig_j
        self.terminus_j = terminus_j
        self.ovl_len: int = ovl_len
    # end def __init__

    def __repr__(self) -> str:
        return '<{}-{}; {}-{}; len={}>'\
            .format(self.contig_i, self.terminus_i,
                self.contig_j, self.terminus_j, self.ovl_len)
    # end def __repr__

    def __eq__(self, other) -> bool:
        if self.contig_i != other.contig_i:
            return False
        # end if
        if self.contig_j != other.contig_j:
            return False
        # end if

        if self.terminus_i != other.terminus_i:
            return False
        # end if
        if self.terminus_j != other.terminus_j:
            return False
        # end if

        if self.ovl_len != other.ovl_len:
            return False
        # end if

        return True
    # end def __eq__

    def __hash__(self) -> int:
        return hash(repr(self))
    # end def __hash__

        return True
    # end def __eq__
# end class Overlap


class OverlapCollection:
    # Class represents collection of `Overlap`s.

    def __init__(self) -> None:
        # Our `_collection` will be a dictionary for the sake of constant access time.
        self._collection: Dict[List[Overlap]] = dict()
    # end def

    def __getitem__(self, key: ContigIndex) -> List[Overlap]:
        # Returns list of overlaps associate with `key` contig.
        try:
            item_to_return = self._collection[key]
        except KeyError:
            # Return tuple just for efficiecy,
            #   though, strictly speaking, it should have been a list.
            item_to_return = tuple()
        # end try
        return item_to_return
    # end def __getitem__

    def __len__(self):
        return len(self._collection.keys())
    # end def __len__

    def add_overlap(self, key: ContigIndex, overlap: Overlap) -> None:
        # Function adds overlap to proper list.
        #
        # :param key: key of contig of interest;
        # :param overlap: `Overlap` instance of the overlap to add;
        try:
            self._collection[key].append(overlap) # just append
        except KeyError:
            # If it's the 1-st overlap, initialize the list.
            self._collection[key] = [overlap]
        # end try
    # end def add_overlap

    def __repr__(self) -> str:
        return str(self._collection)
    # end def __repr__
# end class OverlapCollection


def detect_adjacent_contigs(contig_collection: ContigCollection,
                            mink: int, maxk: int) -> OverlapCollection:
    # Function detects adjacent contigs by comparing their termini.
    #
    # :param contig_collection: instance of ContigCollection returned by
    #   `src.contigs.get_contig_collection` function;
    # :param mink: minimum length of and overlap to be detected;
    # :param maxk: maximum length of and overlap to be detected;

    # Count contigs and save this length in order nom to re-cont it later.
    num_contigs: int = len(contig_collection)

    # Initialize `OverlapCollection` instance
    overlap_collection: OverlapCollection = OverlapCollection()

    # Iterate over contigs and compare it's termini to other termini
    ovl_len: int
    i: ContigIndex
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
        if not ovl_len in (0, contig_collection[i].length):
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
        # We do it in order not to compare pairs of contigs more than one time
        j: ContigIndex
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
                overlap_collection.add_overlap(j, Overlap(j, START, i, RCSTART, ovl_len))
            # end if

            # === Compare i-th end to reverse-complement j-th end ===
            ovl_len = find_overlap_e2s(contig_collection[i].end,
                                       contig_collection[j].rcend,
                                       mink, maxk)
            if ovl_len != 0:
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCEND, ovl_len))
                overlap_collection.add_overlap(j, Overlap(j, END, i, RCEND, ovl_len))
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
