# -*- coding: utf-8 -*-

import pytest

import src.expected_genome_size as egs

from tests.mock_contigs import mock_contigs


@pytest.fixture
def mock_matrix():
    matrix = egs.MatchMatrix()
    matrix[(1, 2)] = 3
    return matrix
# end def mock_matrix

class TestMatchMatrix:

    def test_getitem(self, mock_matrix):
        assert mock_matrix[(1, 2)] == 3
    # end def test_additem_getitem

    def test_getitem_zero(self, mock_matrix):
        assert mock_matrix[(2, 3)] == 0
    # end def test_additem_getitem

    def test_get_overlaps_by_index(self, mock_matrix):
        mock_matrix[(1, 3)] = 4
        assert mock_matrix.get_overlaps_by_index(1) == (((1, 2), 3), ((1, 3), 4))
    # end def test_get_overlaps_by_index
# end class TestMatchMatrix


@pytest.fixture
def mock_matrices(mock_contigs):
    contig_collection, overlap_collection = mock_contigs
    return egs._fill_matrices(overlap_collection)
# end def mock_matrices


class TestExpectedGenomeSize:

    def test_fill_matrices(self, mock_matrices):

        assert repr(mock_matrices[0]) == '{(0, 1): 12, (1, 0): 13, (4, 6): 8}'
        assert repr(mock_matrices[1]) == '{(0, 1): 13, (1, 0): 12, (6, 4): 8}'
        assert repr(mock_matrices[2]) == '{(4, 5): 14, (5, 4): 14}'
        assert repr(mock_matrices[3]) == '{(2, 3): 17, (3, 2): 17}'
    # end def test_fill_matrices

    def test_traceback(self, mock_contigs):
        contig_collection, overlap_collection = mock_contigs
        matrices = egs._fill_matrices(overlap_collection)
        assert egs._trace_back(matrices, contig_collection) == 414
    # end def test_traceback


# end class TestExpectedGenomeSize

