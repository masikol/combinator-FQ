# -*- coding: utf-8 -*-

import pytest
from typing import Tuple

import src.find_overlap as fov


FixtureForFindOvl = Tuple[str, str, int, int]


# === Fixtures for testing `src.find_overlap.find_overlap_s2s` ===

@pytest.fixture
def s2s_ovl_7_7() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 7
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #|||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def s2s_ovl_7_7

@pytest.fixture
def s2s_ovl_7_6() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 6
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #||||||
        'GATCGAAGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def s2s_ovl_7_6

@pytest.fixture
def s2s_ovl_7_8() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 8
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #||||||||
        'GATCGAGTTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def s2s_ovl_7_8

@pytest.fixture
def s2s_ovl_1_20_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 20, overlap length = 7
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #|||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        1,                                     # mink
        20                                     # maxk
    )
# end def s2s_ovl_1_20_7

@pytest.fixture
def s2s_ovl_1_6_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 6, overlap length = 7
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #|||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        1,                                     # mink
        6                                      # maxk
    )
# end def s2s_ovl_1_6_7

@pytest.fixture
def s2s_ovl_8_20_7() -> FixtureForFindOvl:
    # mink = 8, maxk = 20, overlap length = 7
    return (
        'GATCGAGTGTACAGTGAACAATGCTAGGGAGAGCT', # seq1
        #|||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        8,                                     # mink
        20                                     # maxk
    )
# end def s2s_ovl_8_20_7


# === Fixtures for testing `src.find_overlap.find_overlap_e2s` ===

@pytest.fixture
def e2s_ovl_7_7() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 7
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGAG',                             # seq1
        #                            |||||||
                                    'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                                                 # mink
        7                                                                  # maxk
    )
# end def e2s_ovl_7_7

@pytest.fixture
def e2s_ovl_7_6() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 6
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGAG',                             # seq1
        #                            ||||||
                                    'GATCGAAGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                                                 # mink
        7                                                                  # maxk
    )
# end def e2s_ovl_7_6

@pytest.fixture
def e2s_ovl_7_8() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 8
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGAG',                             # seq1
        #                           ||||||||
                                   'GGATCGAGTGAAGAGCCCTAATGTGTAAAATTAAT',  # seq2
        7,                                                                 # mink
        7                                                                  # maxk
    )
# end def e2s_ovl_7_8

@pytest.fixture
def e2s_ovl_1_20_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 20, overlap length = 7
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGAG',                             # seq1
        #                            |||||||
                                    'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        1,                                                                 # mink
        20                                                                 # maxk
    )
# end def e2s_ovl_1_20_7

@pytest.fixture
def e2s_ovl_1_6_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 6, overlap length = 7
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGCT',                             # seq1
        #                            |||||||
                                    'GATCGCTGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        1,                                                                 # mink
        6                                                                  # maxk
    )
# end def e2s_ovl_1_6_7

@pytest.fixture
def e2s_ovl_8_20_7() -> FixtureForFindOvl:
    # mink = 8, maxk = 20, overlap length = 7
    return (
        'GAGAGCTTGTACAGTGAACAATGCTAGGGATCGCT',                             # seq1
        #                            |||||||
                                    'GATCGCTGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        8,                                                                 # mink
        20                                                                 # maxk
    )
# end def e2s_ovl_8_20_7


# === Fixtures for testing `src.find_overlap.find_overlap_e2e` ===

@pytest.fixture
def e2e_ovl_7_7() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 7
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                            |||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAAATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def e2e_ovl_7_7

@pytest.fixture
def e2e_ovl_7_6() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 6
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                             ||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAAGATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def e2e_ovl_7_6

@pytest.fixture
def e2e_ovl_7_8() -> FixtureForFindOvl:
    # Exact k = 7, overlap length = 8
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                           ||||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTAGAATTAAT', # seq2
        7,                                     # mink
        7                                      # maxk
    )
# end def e2e_ovl_7_8

@pytest.fixture
def e2e_ovl_1_20_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 20, overlap length = 7
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                            |||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTACAATTAAT', # seq2
        1,                                     # mink
        20                                     # maxk
    )
# end def e2e_ovl_1_20_7

@pytest.fixture
def e2e_ovl_1_6_7() -> FixtureForFindOvl:
    # mink = 1, maxk = 6, overlap length = 7
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                            |||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTACAATTAAT', # seq2
        1,                                     # mink
        6                                      # maxk
    )
# end def e2e_ovl_1_6_7

@pytest.fixture
def e2e_ovl_8_20_7() -> FixtureForFindOvl:
    # mink = 8, maxk = 20, overlap length = 7
    return (
        'GAGCCCTAGTACAGTGAACAATGCTAGGAATTAAT', # seq1
        #                            |||||||
        'GATCGAGGTGAAGAGCCCTAATGTGTACAATTAAT', # seq2
        8,                                     # mink
        20                                     # maxk
    )
# end def e2e_ovl_8_20_7

# === Testing classes ===

class TestS2S:
    # Test class for `src.find_overlap.find_overlap_s2s`

    def test_s2s_ovl_7_7(self, s2s_ovl_7_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = s2s_ovl_7_7[0]
        seq2: str = s2s_ovl_7_7[1]
        mink: int = s2s_ovl_7_7[2]
        maxk: int = s2s_ovl_7_7[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 7
    # end def test_s2s_ovl_7_7

    def test_s2s_ovl_7_6(self, s2s_ovl_7_6: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = s2s_ovl_7_6[0]
        seq2: str = s2s_ovl_7_6[1]
        mink: int = s2s_ovl_7_6[2]
        maxk: int = s2s_ovl_7_6[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 0
    # end def test_s2s_ovl_7_6

    def test_s2s_ovl_7_8(self, s2s_ovl_7_8: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = s2s_ovl_7_8[0]
        seq2: str = s2s_ovl_7_8[1]
        mink: int = s2s_ovl_7_8[2]
        maxk: int = s2s_ovl_7_8[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 7
    # end def test_s2s_ovl_7_8

    def test_s2s_ovl_1_20_7(self, s2s_ovl_1_20_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = s2s_ovl_1_20_7[0]
        seq2: str = s2s_ovl_1_20_7[1]
        mink: int = s2s_ovl_1_20_7[2]
        maxk: int = s2s_ovl_1_20_7[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 7
    # end def test_s2s_ovl_7_8

    def test_s2s_ovl_1_6_7(self, s2s_ovl_1_6_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = s2s_ovl_1_6_7[0]
        seq2: str = s2s_ovl_1_6_7[1]
        mink: int = s2s_ovl_1_6_7[2]
        maxk: int = s2s_ovl_1_6_7[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 6
    # end def test_s2s_ovl_1_6_7

    def test_s2s_ovl_8_20_7(self, s2s_ovl_8_20_7: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = s2s_ovl_8_20_7[0]
        seq2: str = s2s_ovl_8_20_7[1]
        mink: int = s2s_ovl_8_20_7[2]
        maxk: int = s2s_ovl_8_20_7[3]

        assert fov.find_overlap_s2s(seq1, seq2, mink, maxk) == 0
    # end def test_s2s_ovl_1_6_7
# end class TestS2S


class TestE2S:
    # Test class for `src.find_overlap.find_overlap_e2s`

    def test_e2s_ovl_7_7(self, e2s_ovl_7_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2s_ovl_7_7[0]
        seq2: str = e2s_ovl_7_7[1]
        mink: int = e2s_ovl_7_7[2]
        maxk: int = e2s_ovl_7_7[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 7
    # end def test_e2s_ovl_7_7

    def test_e2s_ovl_7_6(self, e2s_ovl_7_6: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = e2s_ovl_7_6[0]
        seq2: str = e2s_ovl_7_6[1]
        mink: int = e2s_ovl_7_6[2]
        maxk: int = e2s_ovl_7_6[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 0
    # end def test_e2s_ovl_7_6

    def test_e2s_ovl_7_8(self, e2s_ovl_7_8: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = e2s_ovl_7_8[0]
        seq2: str = e2s_ovl_7_8[1]
        mink: int = e2s_ovl_7_8[2]
        maxk: int = e2s_ovl_7_8[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 0
    # end def test_e2s_ovl_7_8

    def test_e2s_ovl_1_20_7(self, e2s_ovl_1_20_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2s_ovl_1_20_7[0]
        seq2: str = e2s_ovl_1_20_7[1]
        mink: int = e2s_ovl_1_20_7[2]
        maxk: int = e2s_ovl_1_20_7[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 7
    # end def test_e2s_ovl_1_20_7

    def test_e2s_ovl_1_6_7(self, e2s_ovl_1_6_7: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = e2s_ovl_1_6_7[0]
        seq2: str = e2s_ovl_1_6_7[1]
        mink: int = e2s_ovl_1_6_7[2]
        maxk: int = e2s_ovl_1_6_7[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 0
    # end def test_e2s_ovl_1_6_7

    def test_e2s_ovl_8_20_7(self, e2s_ovl_8_20_7: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = e2s_ovl_8_20_7[0]
        seq2: str = e2s_ovl_8_20_7[1]
        mink: int = e2s_ovl_8_20_7[2]
        maxk: int = e2s_ovl_8_20_7[3]

        assert fov.find_overlap_e2s(seq1, seq2, mink, maxk) == 0
    # end def test_e2s_ovl_8_20_7
# end class TestE2S


class TestE2E:
    # Test class for `src.find_overlap.find_overlap_e2e`

    def test_e2e_ovl_7_7(self, e2e_ovl_7_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2e_ovl_7_7[0]
        seq2: str = e2e_ovl_7_7[1]
        mink: int = e2e_ovl_7_7[2]
        maxk: int = e2e_ovl_7_7[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 7
    # end def test_e2e_ovl_7_7

    def test_e2e_ovl_7_6(self, e2e_ovl_7_6: FixtureForFindOvl):
        # Should not find any overlap
        seq1: str = e2e_ovl_7_6[0]
        seq2: str = e2e_ovl_7_6[1]
        mink: int = e2e_ovl_7_6[2]
        maxk: int = e2e_ovl_7_6[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 0
    # end def test_e2e_ovl_7_6

    def test_e2e_ovl_7_8(self, e2e_ovl_7_8: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2e_ovl_7_8[0]
        seq2: str = e2e_ovl_7_8[1]
        mink: int = e2e_ovl_7_8[2]
        maxk: int = e2e_ovl_7_8[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 7
    # end def test_e2e_ovl_7_8

    def test_e2e_ovl_1_20_7(self, e2e_ovl_1_20_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2e_ovl_1_20_7[0]
        seq2: str = e2e_ovl_1_20_7[1]
        mink: int = e2e_ovl_1_20_7[2]
        maxk: int = e2e_ovl_1_20_7[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 7
    # end def test_e2e_ovl_1_20_7

    def test_e2e_ovl_1_6_7(self, e2e_ovl_1_6_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2e_ovl_1_6_7[0]
        seq2: str = e2e_ovl_1_6_7[1]
        mink: int = e2e_ovl_1_6_7[2]
        maxk: int = e2e_ovl_1_6_7[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 6
    # end def test_e2e_ovl_1_6_7

    def test_e2e_ovl_8_20_7(self, e2e_ovl_8_20_7: FixtureForFindOvl):
        # Should find overlap of 7 bp length
        seq1: str = e2e_ovl_8_20_7[0]
        seq2: str = e2e_ovl_8_20_7[1]
        mink: int = e2e_ovl_8_20_7[2]
        maxk: int = e2e_ovl_8_20_7[3]

        assert fov.find_overlap_e2e(seq1, seq2, mink, maxk) == 0
    # end def test_e2e_ovl_8_20_7
# end class TestE2E
