# -*- coding: utf-8 -*-

import pytest

import src.expected_genome_size as egs

from tests.mock_contigs import mock_contigs


@pytest.fixture
def mock_martix():
    matrix = egs.MatchMatrix()
    matrix[(1, 2)] = 3
    return matrix
# end def mock_martix

class TestMatchMatrix:

    def test_getitem(self, mock_martix):
        assert mock_martix[(1, 2)] == 3
    # end def test_additem_getitem

    def test_getitem_zero(self, mock_martix):
        assert mock_martix[(2, 3)] == 0
    # end def test_additem_getitem
# end class TestMatchMatrix


@pytest.fixture
def mock_matrices(mock_contigs):
    contig_collection, overlap_collection = mock_contigs
    return egs._fill_matrices(overlap_collection)
# end def mock_matrices


class TestExpectedGenomeSize:

    def test_matrix(self, mock_matrices):

        assert repr(mock_matrices[0]) == '{(0, 1): 12, (1, 0): 13, (4, 6): 8}'
        assert repr(mock_matrices[1]) == '{(0, 1): 13, (1, 0): 12, (6, 4): 8}'
        assert repr(mock_matrices[2]) == '{(4, 5): 14, (5, 4): 14}'
        assert repr(mock_matrices[3]) == '{(2, 3): 17, (3, 2): 17}'
    # end def test_matrix


# end class TestExpectedGenomeSize

