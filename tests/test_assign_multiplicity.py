# -*- coding: utf-8 -*-

import os
import pytest
from typing import MutableSequence, Sequence, Tuple

import src.contigs as cnt
import src.overlaps as ovl
import src.assign_multiplicity as amu
from src.overlaps import START, RCSTART, END, RCEND


# === Fixtures for function `src.assign_multiplicity._calc_multiplty_by_coverage` ===
@pytest.fixture
def coverages_equal() -> Tuple[float, float]:
    return 19.89, 19.89
# end def coverages_equal

@pytest.fixture
def coverages_x3_dot_45_repeat() -> Tuple[float, float]:
    first_cov: float = 19.89
    return first_cov * 3.45, 19.89
# end def coverages_x3_dot_45_repeat

@pytest.fixture
def coverages_x0_dot_45_repeat() -> Tuple[float, float]:
    first_cov: float = 19.89
    return first_cov * 0.45, 19.89
# end def coverages_x0_dot_45_repeat


# === Fixtures for function `src.assign_multiplicity._calc_multiplty_by_coverage` ===
@pytest.fixture
def overlaps_empty() -> MutableSequence[ovl.Overlap]:
    return []
# end def overlaps_empty

@pytest.fixture
def overlaps_start1_end0() -> MutableSequence[ovl.Overlap]:
    return [
        ovl.Overlap(1, START, 2, END, 15)
    ]
# end def overlaps_start1_end0

@pytest.fixture
def overlaps_start1_end1() -> MutableSequence[ovl.Overlap]:
    return [
        ovl.Overlap(1, START, 2, END, 15),
        ovl.Overlap(1, END, 2, START, 15)
    ]
# end def overlaps_start1_end1

@pytest.fixture
def overlaps_start1_end2() -> MutableSequence[ovl.Overlap]:
    return [
        ovl.Overlap(1, START, 2, END, 15),
        ovl.Overlap(1, END, 2, START, 15),
        ovl.Overlap(1, END, 3, START, 15)
    ]
# end def overlaps_start1_end2

@pytest.fixture
def overlaps_start2_end2() -> MutableSequence[ovl.Overlap]:
    return [
        ovl.Overlap(1, START, 2, END, 15),
        ovl.Overlap(1, START, 4, RCSTART, 15),
        ovl.Overlap(1, END, 2, START, 15),
        ovl.Overlap(1, END, 3, START, 15)
    ]
# end def overlaps_start2_end2

@pytest.fixture
def overlaps_start2_end4() -> MutableSequence[ovl.Overlap]:
    return [
        ovl.Overlap(1, START, 2, END, 15),
        ovl.Overlap(1, START, 4, RCSTART, 15),
        ovl.Overlap(1, END, 2, START, 15),
        ovl.Overlap(1, END, 3, RCEND, 15),
        ovl.Overlap(1, END, 5, RCEND, 15),
        ovl.Overlap(1, END, 7, START, 15)
    ]
# end def overlaps_start2_end4


# === Fixtures for function `src.assign_multiplicity.assign_multiplty` ===

@pytest.fixture
def contig_collection_spades_0() -> cnt.ContigCollection:

    mink: int = 16
    maxk: int = 25

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_spades_0

@pytest.fixture
def contig_collection_spades_1() -> cnt.ContigCollection:

    mink: int = 8
    maxk: int = 17

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_1.fasta.gz'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_spades_1

@pytest.fixture
def contig_collection_a5_0() -> cnt.ContigCollection:

    mink: int = 16
    maxk: int = 25

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_a5_0

@pytest.fixture
def contig_collection_mix_0() -> cnt.ContigCollection:

    mink: int = 16
    maxk: int = 25

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_mix_0.fasta'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_mix_0

@pytest.fixture
def contig_collection_spades_no_multplty_0() -> cnt.ContigCollection:

    mink: int = 16
    maxk: int = 25

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_no_multplty_0.fasta'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_spades_no_multplty_0

@pytest.fixture
def contig_collection_a5_repeat() -> cnt.ContigCollection:

    mink: int = 15
    maxk: int = 15

    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_a5_repeat.fasta'),
        maxk
    )

    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def contig_collection_a5_repeat


class TestCalcMultpltyByCov:
    # Class for testing function `src.assign_multiplicity._calc_multiplty_by_coverage`

    def test_calc_multplty_by_cov_equal(self, coverages_equal):
        # Test how the function calculates multiplicity for equal coverages.

        curr_cov: float = coverages_equal[0]
        first_cov: float = coverages_equal[1]

        expected: float = 1.0

        obtained: float = amu._calc_multiplty_by_coverage(curr_cov, first_cov)

        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_cov_equal

    def test_calc_multplty_by_cov_x3_dot_45_repeat(self, coverages_x3_dot_45_repeat):
        # Test how the function calculates multiplicity for contig
        #   of coverage 3.45 times greater than coverage of the first contig.

        curr_cov: float = coverages_x3_dot_45_repeat[0]
        first_cov: float = coverages_x3_dot_45_repeat[1]

        expected: float = 3.45

        obtained: float = amu._calc_multiplty_by_coverage(curr_cov, first_cov)

        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_cov_x3_dot_45_repeat

    def test_calc_multplty_by_cov_x0_dot_45_repeat(self, coverages_x0_dot_45_repeat):
        # Test how the function calculates multiplicity for contig
        #   of coverage 0.45 times "greater" than coverage of the first contig.

        curr_cov: float = coverages_x0_dot_45_repeat[0]
        first_cov: float = coverages_x0_dot_45_repeat[1]

        expected: float = 1.0

        obtained: float = amu._calc_multiplty_by_coverage(curr_cov, first_cov)

        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_cov_x0_dot_45_repeat
# end class TestCalcMultpltyByCov


class TestCalcMultpltyByOvl:
    # Class for testing function `src.assign_multiplicity._calc_multiplty_by_overlaps`

    def test_calc_multplty_by_ovl_empty(self, overlaps_empty):
        # Tests how the function assigns mulpitlicity if
        #   contig has no overlaps.
        expected: float = 1.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_empty)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_empty

    def test_calc_multplty_by_ovl_start1_end0(self, overlaps_start1_end0):
        # Tests how the function assigns mulpitlicity if
        #   contig has single overlap at one terminus.
        expected: float = 1.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_start1_end0)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_start1_end0

    def test_calc_multplty_by_ovl_start1_end1(self, overlaps_start1_end1):
        # Tests how the function assigns mulpitlicity if
        #   contig has one overlap at both termini.
        expected: float = 1.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_start1_end1)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_start1_end1

    def test_calc_multplty_by_ovl_start1_end2(self, overlaps_start1_end2):
        # Tests how the function assigns mulpitlicity if
        #   contig has one overlap at one terminus
        #   and two overlaps -- at another.
        expected: float = 1.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_start1_end2)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_start1_end2

    def test_calc_multplty_by_ovl_start2_end2(self, overlaps_start2_end2):
        # Tests how the function assigns mulpitlicity if
        #   contig has two overlaps at both termini.
        expected: float = 2.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_start2_end2)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_start2_end2

    def test_calc_multplty_by_ovl_start2_end4(self, overlaps_start2_end4):
        # Tests how the function assigns mulpitlicity if
        #   contig has two overlaps at one terminus
        #   and 4 overlaps -- at another.
        expected: float = 2.0
        obtained: float = amu._calc_multiplty_by_overlaps(overlaps_start2_end4)
        assert abs(obtained - expected) < 1e-1
    # end def test_calc_multplty_by_ovl_start2_end4

# end class TestCalcMultpltyByOvl


class TestAssignMultiplty:
    # Class for testing function `src.assign_multiplicity.assign_multiplty`

    def test_assign_multiplty_spades_0(self, contig_collection_spades_0):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_spades_0.fasta

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_spades_0

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = tuple([
            1.0, # NODE_1
            1.0, # NODE_2
            1.8, # NODE_3
            1.2  # NODE_4
        ])

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_spades_0

    def test_assign_multiplty_spades_1(self, contig_collection_spades_1):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_spades_1.fasta.gz`

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_spades_1

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = tuple([1.0] * 7)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_spades_1

    def test_assign_multiplty_a5_0(self, contig_collection_a5_0):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_a5_0.fasta`

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_a5_0

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = tuple([1.0] * 4)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_a5_0

    def test_assign_multiplty_mix_0(self, contig_collection_mix_0):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_mix_0.fasta`

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_mix_0

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = tuple([1.0] * 4)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_mix_0

    def test_assign_multiplty_spades_no_multplty_0(self, contig_collection_spades_no_multplty_0):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_spades_no_multplty_0.fasta`

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_spades_no_multplty_0

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = tuple([1.0] * 4)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_spades_no_multplty_0

    def test_assign_multiplty_a5_repeat(self, contig_collection_a5_repeat):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_spades_no_multplty_0.fasta`

        # "Rename" variable
        contig_collection: cnt.ContigCollection
        overlap_collection: ovl.OverlapCollection
        contig_collection, overlap_collection = contig_collection_a5_repeat

        # Assign
        amu.assign_multiplty(contig_collection, overlap_collection)

        expected_multplts: Sequence[float] = (1.0, 1.0, 2.0, 1.0, 1.0)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_a5_repeat
# end class TestAssignMultiplty
