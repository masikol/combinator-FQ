# -*- coding: utf-8 -*-

import os
import pytest
from typing import List, Generator

import src.contigs as cnt
import src.overlaps as ovl
from src.overlaps import START, RCSTART, END, RCEND

from tests.mock_contigs import mock_contigs_spades_0


# === Fixtures for testing function `src.overlaps.detect_adjacent_contigs` ===

@pytest.fixture
def contig_collection_spades_0() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
        25 # maxk
    )

    return contig_collection, 16, 25
# end def contig_collection_spades_0


# === Test classes ===


class TestOverlap:
    # Class for testing class `src.overlaps.Overlap`

    ovl1: ovl.Overlap = ovl.Overlap(0, END, 1, START, 33)

    @staticmethod
    def _ovl2_generator() -> Generator[ovl.Overlap, None, None]:
        yield ovl.Overlap(0, END, 1, START, 33) # equals to ovl1
        yield ovl.Overlap(1, END, 1, START, 33)
        #                 ^
        yield ovl.Overlap(1, RCSTART, 1, START, 33)
        #                    ^
        yield ovl.Overlap(1, END, 2, START, 33)
        #                         ^
        yield ovl.Overlap(1, END, 1, RCEND, 33)
        #                            ^
        yield ovl.Overlap(1, END, 1, START, 21)
        #                                   ^
    # end def _ovl2_generator

    def test_eq(self):
        # Tests method __eq__

        ovl2_iterator = self._ovl2_generator()
        ovl2: ovl.Overlap = next(ovl2_iterator)
        assert self.ovl1 == ovl2 # first element should be equal

        for ovl2 in ovl2_iterator:
            assert self.ovl1 != ovl2
        # end for
    # end def test_eq

    def test_hash(self):
        # Tests method __hash__

        ovl2_iterator = self._ovl2_generator()
        ovl2: ovl.Overlap = next(ovl2_iterator)
        assert hash(self.ovl1) == hash(ovl2) # first element should be equal

        for ovl2 in ovl2_iterator:
            assert hash(self.ovl1) != hash(ovl2)
        # end for
    # end def test_hash
# end class TestOverlap


class TestOverlapCollection:
    # Class for testing class `src.overlaps.OverlapCollection`

    test_instance: ovl.OverlapCollection = ovl.OverlapCollection()
    first_index: cnt.ContigIndex = 0
    second_index: cnt.ContigIndex = 1
    nonextant_index: cnt.ContigIndex = 3

    def test_add_overlap(self):
        # Function for testing method `OverlapCollection.add_overlap`
        self.test_instance.add_overlap(self.first_index, ovl.Overlap(self.first_index, START, 1, END, 21))
        try:
            assert len(self.test_instance[self.first_index]) == 1
        except KeyError:
            pytest.fatal('Error: method `OverlapCollection.add_overlap` does not work.`')
        # end try

        self.test_instance.add_overlap(self.second_index, ovl.Overlap(self.second_index, START, 4, END, 55))
        try:
            assert len(self.test_instance) == 2
        except KeyError:
            pytest.fatal('Error: method `OverlapCollection.add_overlap` does not work.`')
        # end try

        self.test_instance.add_overlap(self.first_index, ovl.Overlap(self.first_index, START, 3, END, 33))
        try:
            assert len(self.test_instance[self.first_index]) == 2
        except KeyError:
            pytest.fatal('Error: method `OverlapCollection.add_overlap` does not work.`')
        # end try
    # end def test_add_overlap

    def test_getitem(self):
        # Function for testing method `OverlapCollection.__getitem__`
        try:
            assert len(self.test_instance[self.first_index]) == 2
            assert len(self.test_instance[self.nonextant_index]) == 0
        except KeyError:
            pytest.fatal('Error: method `OverlapCollection.__getitem__` does not work.`')
        # end try
    # end def test_getitem
# end class TestOverlapCollection


class TestDetectAdjacentContigs:
    # Class for testing function `src.overlaps.detect_adjacent_contigs`

    def test_detect_adjacent_contigs(self, contig_collection_spades_0):
        contig_collection, mink, maxk = contig_collection_spades_0
        overlap_collection: cnt.OverlapCollection = ovl.detect_adjacent_contigs(
            contig_collection,
            mink,
            maxk
        )

        node_1: ContigIndex = 0
        node_2: ContigIndex = 1
        node_3: ContigIndex = 2
        node_4: ContigIndex = 3

        expected: List[ovl.Overlap]

        expected = {
            ovl.Overlap(node_1, START, node_2, END, 18),
            ovl.Overlap(node_1, END, node_2, START, 16)
        }
        assert set(overlap_collection[node_1]) == expected

        expected = {
            ovl.Overlap(node_2, START, node_1, END, 16),
            ovl.Overlap(node_2, END, node_1, START, 18)
        }
        assert set(overlap_collection[node_2]) == expected

        expected = {
            ovl.Overlap(node_3, START, node_3, END, 21),
            ovl.Overlap(node_3, END, node_3, START, 21)
        }
        assert set(overlap_collection[node_3]) == expected

        expected = set()
        assert set(overlap_collection[node_4]) == expected

    # end def test_detect_adjacent_contigs
# end class TestDetectAdjacentContigs