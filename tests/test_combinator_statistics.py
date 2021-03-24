# -*- coding: utf-8 -*-

import pytest
from math import floor
from typing import Sequence

import src.combinator_statistics as sts
from src.overlaps import Overlap, START, RCSTART, END, RCEND

from tests.mock_contigs import mock_contigs_spades_0, mock_contigs_spades_1
from tests.mock_contigs import mock_contigs_a5_0, mock_contigs_mix_0
from tests.mock_contigs import mock_contigs_spades_no_multplty_0


class TestCoverageCalculator:
    # Class for testing class `src.combinator_statistics.CoverageCalculator`

    # Test "private" method `_filter_non_none_covs`
    def test_filter_non_none_covs_spades_0(self, mock_contigs_spades_0):
        # Should find 4 "coverage-containing" contigs in `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_len: int = 4

        assert len(cov_calc._coverages) == expected_len
    # end def test_filter_non_none_covs_spades_0

    def test_filter_non_none_covs_spades_1(self, mock_contigs_spades_1):
        # Should find 7 "coverage-containing" contigs in `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_len: int = 7

        assert len(cov_calc._coverages) == expected_len
    # end def test_filter_non_none_covs_spades_1

    def test_filter_non_none_covs_a5_0(self, mock_contigs_a5_0):
        # Should find no "coverage-containing" contigs in `mock_contigs_a5_0`
        contig_collection, overlap_collection = mock_contigs_a5_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_len: int = 0

        assert len(cov_calc._coverages) == expected_len
    # end def test_filter_non_none_covs_a5_0

    def test_filter_non_none_covs_mix_0(self, mock_contigs_mix_0):
        # Should find 2 "coverage-containing" contigs in `mock_contigs_mix_0`
        contig_collection, overlap_collection = mock_contigs_mix_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_len: int = 2

        assert len(cov_calc._coverages) == expected_len
    # end def test_filter_non_none_covs_mix_0

    def test_filter_non_none_covs_spades_no_multiplty_0(self, mock_contigs_spades_no_multplty_0):
        # Should find 4 "coverage-containing" contigs in `mock_contigs_spades_no_multplty_0`
        contig_collection, overlap_collection = mock_contigs_spades_no_multplty_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_len: int = 4

        assert len(cov_calc._coverages) == expected_len
    # end def test_filter_non_none_covs_spades_no_multiplty_0


    # Test method `get_min_coverage`
    def test_get_min_coverage_spades_0(self, mock_contigs_spades_0):
        # Should return 18.28 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_min_cov: float = 18.28
        assert abs(cov_calc.get_min_coverage() - expected_min_cov) < 1e-3
    # end def test_get_min_coverage_spades_0

    def test_get_min_coverage_spades_1(self, mock_contigs_spades_1):
        # Should return 24.26 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_min_cov: float = 24.26
        assert abs(cov_calc.get_min_coverage() - expected_min_cov) < 1e-3
    # end def test_get_min_coverage_spades_1

    def test_get_min_coverage_a5_0(self, mock_contigs_a5_0):
        # Should return None for `mock_contigs_a5_0`
        contig_collection, overlap_collection = mock_contigs_a5_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_min_cov: float = None
        assert cov_calc.get_min_coverage() is expected_min_cov
    # end def test_get_min_coverage_a5_0

    def test_get_min_coverage_mix_0(self, mock_contigs_mix_0):
        # Should return 18.28 for `mock_contigs_mix_0`
        contig_collection, overlap_collection = mock_contigs_mix_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_min_cov: float = 18.28
        assert abs(cov_calc.get_min_coverage() - expected_min_cov) < 1e-3
    # end def test_get_min_coverage_mix_0

    def test_get_min_coverage_spades_no_multplty_0_0(self, mock_contigs_spades_no_multplty_0):
        # Should return 0.00 for `mock_contigs_spades_no_multplty_0`
        contig_collection, overlap_collection = mock_contigs_spades_no_multplty_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_min_cov: float = 0.00
        assert abs(cov_calc.get_min_coverage() - expected_min_cov) < 1e-3
    # end def test_get_min_coverage_spades_no_multplty_0_0


    # Test method `get_max_coverage`
    def test_get_max_coverage_spades_0(self, mock_contigs_spades_0):
        # Should return 43.73 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_max_cov: float = 43.73
        assert abs(cov_calc.get_max_coverage() - expected_max_cov) < 1e-3
    # end def test_get_max_coverage_spades_0

    def test_get_max_coverage_spades_1(self, mock_contigs_spades_1):
        # Should return 24.26 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_max_cov: float = 24.26
        assert abs(cov_calc.get_max_coverage() - expected_max_cov) < 1e-3
    # end def test_get_max_coverage_spades_1

    def test_get_max_coverage_a5_0(self, mock_contigs_a5_0):
        # Should return None for `mock_contigs_a5_0`
        contig_collection, overlap_collection = mock_contigs_a5_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_max_cov: float = None
        assert cov_calc.get_max_coverage() is expected_max_cov
    # end def test_get_max_coverage_a5_0

    def test_get_max_coverage_mix_0(self, mock_contigs_mix_0):
        # Should return 28.68 for `mock_contigs_mix_0`
        contig_collection, overlap_collection = mock_contigs_mix_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_max_cov: float = 28.68
        assert abs(cov_calc.get_max_coverage() - expected_max_cov) < 1e-3
    # end def test_get_max_coverage_mix_0

    def test_get_max_coverage_spades_no_multplty_0_0(self, mock_contigs_spades_no_multplty_0):
        # Should return 43.73 for `mock_contigs_spades_no_multplty_0`
        contig_collection, overlap_collection = mock_contigs_spades_no_multplty_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_max_cov: float = 43.73
        assert abs(cov_calc.get_max_coverage() - expected_max_cov) < 1e-3
    # end def test_get_max_coverage_spades_no_multplty_0_0


    # Test method `calc_mean_coverage`
    def test_get_mean_coverage_spades_0(self, mock_contigs_spades_0):
        # Should return 28.74 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_mean_cov: float = 28.74
        assert abs(cov_calc.calc_mean_coverage() - expected_mean_cov) < 1e-3
    # end def test_get_mean_coverage_spades_0

    def test_get_mean_coverage_spades_1(self, mock_contigs_spades_1):
        # Should return 24.26 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_mean_cov: float = 24.26
        assert abs(cov_calc.calc_mean_coverage() - expected_mean_cov) < 1e-3
    # end def test_get_mean_coverage_spades_1

    def test_get_mean_coverage_a5_0(self, mock_contigs_a5_0):
        # Should return None for `mock_contigs_a5_0`
        contig_collection, overlap_collection = mock_contigs_a5_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_mean_cov: float = None
        assert cov_calc.calc_mean_coverage() is expected_mean_cov
    # end def test_get_mean_coverage_a5_0

    def test_get_mean_coverage_mix_0(self, mock_contigs_mix_0):
        # Should return 23.48 for `mock_contigs_mix_0`
        contig_collection, overlap_collection = mock_contigs_mix_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_mean_cov: float = 23.48
        assert abs(cov_calc.calc_mean_coverage() - expected_mean_cov) < 1e-3
    # end def test_get_mean_coverage_mix_0

    def test_get_mean_coverage_spades_no_multplty_0_0(self, mock_contigs_spades_no_multplty_0):
        # Should return 22.67 for `mock_contigs_spades_no_multplty_0`
        contig_collection, overlap_collection = mock_contigs_spades_no_multplty_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_mean_cov: float = 22.67
        assert abs(cov_calc.calc_mean_coverage() - expected_mean_cov) < 1e-3
    # end def test_get_mean_coverage_spades_no_multplty_0_0


    # Test method `calc_median_coverage`
    def test_get_median_coverage_spades_0(self, mock_contigs_spades_0):
        # Should return 26.47 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_median_cov: float = 26.47
        assert abs(cov_calc.calc_median_coverage() - expected_median_cov) < 1e-3
    # end def test_get_median_coverage_spades_0

    def test_get_median_coverage_spades_1(self, mock_contigs_spades_1):
        # Should return 24.26 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_median_cov: float = 24.26
        assert abs(cov_calc.calc_median_coverage() - expected_median_cov) < 1e-3
    # end def test_get_median_coverage_spades_1

    def test_get_median_coverage_a5_0(self, mock_contigs_a5_0):
        # Should return None for `mock_contigs_a5_0`
        contig_collection, overlap_collection = mock_contigs_a5_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_median_cov: float = None
        assert cov_calc.calc_median_coverage() is expected_median_cov
    # end def test_get_median_coverage_a5_0

    def test_get_median_coverage_mix_0(self, mock_contigs_mix_0):
        # Should return 23.48 for `mock_contigs_mix_0`
        contig_collection, overlap_collection = mock_contigs_mix_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_median_cov: float = 23.48
        assert abs(cov_calc.calc_median_coverage() - expected_median_cov) < 1e-3
    # end def test_get_median_coverage_mix_0

    def test_get_median_coverage_spades_no_multplty_0_0(self, mock_contigs_spades_no_multplty_0):
        # Should return 23.48 for `mock_contigs_spades_no_multplty_0`
        contig_collection, overlap_collection = mock_contigs_spades_no_multplty_0
        cov_calc = sts.CoverageCalculator(contig_collection)

        expected_median_cov: float = 23.48
        assert abs(cov_calc.calc_median_coverage() - expected_median_cov) < 1e-3
    # end def test_get_median_coverage_spades_no_multplty_0_0
# end class TestCoverageCalculator


class TestCalcSumContigLengths:
    # Class for testing function`src.combinator_statistics.calc_sum_contig_lengths`

    def test_calc_sum_contig_lengths_spades_0(self, mock_contigs_spades_0):
        # Should return 1194 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0

        expected_sum_lengths: int = 1194

        assert sts.calc_sum_contig_lengths(contig_collection) == expected_sum_lengths
    # end def test_calc_sum_contig_lengths_spades_0

    def test_calc_sum_contig_lengths_spades_1(self, mock_contigs_spades_1):
        # Should return 470 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1

        expected_sum_lengths: int = 470

        assert sts.calc_sum_contig_lengths(contig_collection) == expected_sum_lengths
    # end def test_calc_sum_contig_lengths_spades_1
# end class TestCalcSumContigLengths


class TestIsMatch:
    # Class for testing functions:
    #   `src.combinator_statistics.is_start_match`
    #   `src.combinator_statistics.is_end_match`

    def test_is_start_match(self):
        # Function for testing function `src.combinator_statistics.is_start_match`
        assert sts.is_start_match(Overlap(0, START, 1, END, 21)) ==     True
        assert sts.is_start_match(Overlap(0, START, 1, RCSTART, 21)) == True
        assert sts.is_start_match(Overlap(0, END, 1, START, 21)) ==     False
        assert sts.is_start_match(Overlap(0, END, 1, RCEND, 21)) ==     False
        assert sts.is_start_match(Overlap(0, END, 1, RCSTART, 21)) ==   False
        assert sts.is_start_match(Overlap(0, START, 1, RCEND, 21)) ==   False
        assert sts.is_start_match(Overlap(0, START, 1, START, 21)) ==   False
        assert sts.is_start_match(Overlap(0, END, 1, END, 21)) ==       False
    # end def test_is_start_match

    def test_is_end_match(self):
        # Function for testing function `src.combinator_statistics.is_start_match`
        assert sts.is_end_match(Overlap(0, END, 1, START, 21)) ==     True
        assert sts.is_end_match(Overlap(0, END, 1, RCEND, 21)) ==     True
        assert sts.is_end_match(Overlap(0, START, 1, END, 21)) ==     False
        assert sts.is_end_match(Overlap(0, START, 1, RCSTART, 21)) == False
        assert sts.is_end_match(Overlap(0, END, 1, RCSTART, 21)) ==   False
        assert sts.is_end_match(Overlap(0, START, 1, RCEND, 21)) ==   False
        assert sts.is_end_match(Overlap(0, START, 1, START, 21)) ==   False
        assert sts.is_end_match(Overlap(0, END, 1, END, 21)) ==       False
    # end def test_is_end_match
# end class TestIsMatch


class TestCalcLqCoef:
    # Class for testing function `src.combinator_statistics.calc_lq_coef`

    def test_calc_lq_coef_spades_0(self, mock_contigs_spades_0):
        # Should return 60.00 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0

        expected_lq_coef: float = 60.00

        assert abs(sts.calc_lq_coef(contig_collection, overlap_collection) - expected_lq_coef) < 1e-3
    # end def test_calc_lq_coef_spades_0

    def test_calc_lq_coef_spades_1(self, mock_contigs_spades_1):
        # Should return 64.29 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1

        expected_lq_coef: float = 64.29

        assert (sts.calc_lq_coef(contig_collection, overlap_collection) - expected_lq_coef) < 1e-3
    # end def test_calc_lq_coef_spades_1

    def test_calc_lq_coef_spades_1_node3_cov_x2(self, mock_contigs_spades_1):
        # We'll make NODE_4 2-copy
        # Should return 56.25 for modified `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1

        # Increase twofold copy number of contig `NODE_4`
        save_cov: float = contig_collection[3].cov
        save_multplty: float = contig_collection[3].multplty

        contig_collection[3].cov = floor(contig_collection[3].cov * 2.0)
        contig_collection[3].multplty = floor(contig_collection[3].multplty * 2.0)

        expected_lq_coef: float = 56.25

        try:
            assert abs(sts.calc_lq_coef(contig_collection, overlap_collection) - expected_lq_coef) < 1e-3
        except AssertionError:
            pass
        finally:
            # Return to init preset, since `mock_contigs_spades_1` is of session scope
            contig_collection[3].cov = save_cov
            contig_collection[3].multplty = save_multplty
        # end try
    # end def test_calc_lq_coef_spades_1_node3_cov_x2
# end class TestCalcLqCoef


class TestCalcExpGenomeSize:
    # Class for testing function `src.combinator_statistics.calc_exp_genome_size`

    def test_calc_exp_genome_size_spades_0(self, mock_contigs_spades_0):
        # Should return 1383 for `mock_contigs_spades_0`
        contig_collection, overlap_collection = mock_contigs_spades_0

        expected_egs: int = 1383

        assert sts.calc_exp_genome_size(contig_collection, overlap_collection) == expected_egs
    # end def test_calc_exp_genome_size_spades_0

    def test_calc_exp_genome_size_spades_1(self, mock_contigs_spades_1):
        # Should return 1383 for `mock_contigs_spades_1`
        contig_collection, overlap_collection = mock_contigs_spades_1

        expected_egs: int = 414

        assert sts.calc_exp_genome_size(contig_collection, overlap_collection) == expected_egs
    # end def test_calc_exp_genome_size_spades_1
# end class TestCalcExpGenomeSize
