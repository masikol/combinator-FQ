# -*- coding: utf-8 -*-

import pytest

import src.combinator_statistics as sts

from tests.mock_contigs import mock_contigs


class TestExpectedGenomeSize:

    def test_calc_exp_genome_size(self, mock_contigs):
        contig_collection, overlap_collection = mock_contigs
        assert sts.calc_exp_genome_size(contig_collection, overlap_collection) == 414
    # end def test_traceback
# end class TestExpectedGenomeSize

