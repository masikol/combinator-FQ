# -*- coding: utf-8 -*-

import os
import pytest
from typing import Tuple, Sequence, Callable, Any, Generator

import src.contigs as cnt


# === Fixtures for function `src.contigs._rc` ===
# And for `src.contigs._calc_gc_content`

@pytest.fixture
def some_sequence() -> str:
    return 'ACCATAATGAATCTGGCT'
# end def some_sequence

@pytest.fixture
def degenerate_sequence() -> str:
    return 'ACRASAATGAATCDGGNT'
# end def degenerate_sequence


# === Fixtures for function `src.contigs._validate_fasta` ===

@pytest.fixture
def valid_fasta() -> Tuple[str, str]:
    # Valid fasta record
    return tuple([
        '>some_sequence',
        'ACTACCAACATGAAGAAACAATTTATAAT'
    ])
# end def valid_fasta

@pytest.fixture
def invalid_fasta_header() -> Tuple[str, str]:
    # Fasta record with invalid header
    return tuple([
        'some_sequence',
        'ACTACCAACATGAAGAAACAATTTATAAT'
    ])
# end def invalid_fasta_header

@pytest.fixture
def invalid_empty_seq() -> Tuple[str, str]:
    # Fasta record with invalid sequence
    return tuple([
        '>some_sequence',
        ''
    ])
# end def invalid_empty_seq

@pytest.fixture
def invalid_fasta_seq() -> Tuple[str, str]:
    # Fasta record with invalid sequence
    return tuple([
        '>some_sequence',
        'ACTA$CAACATGAggtaAAACAAXTTATAAT'
    ])
# end def invalid_fasta_seq


# === Fixtures for generator `src.contigs._fasta_generator` ===

@pytest.fixture
def plain_text_fasta() -> str:
    return os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
# end def plain_text_fasta

@pytest.fixture
def gzipped_fasta() -> str:
    return os.path.join('tests', 'data', 'test_contigs_spades_1.fasta.gz')
# end def gzipped_fasta


# === Fixtures for function `src.contigs._is_spades_header` ===
# And for `src.contigs._format_contig_name`
# And for `src.contigs._parse_coverage`

@pytest.fixture
def valid_spades_header() -> str:
    # Returns valid SPAdes header
    return 'NODE_1_length_938688_cov_24.261358'
# end def valid_spades_header

@pytest.fixture
def cryptic_invalid_spades_header() -> str:
    # Returns cryptic invalid SPAdes header
    return 'NODE_1_length_938688_cov_24.26.1,358'
# end def cryptic_invalid_spades_header

@pytest.fixture
def valid_a5_header() -> str:
    # Returns valid a5 header
    return 'scaffold_24'
# end def valid_a5_header


# === Fixtures for function `src.contigs.get_contig_collection` ===

@pytest.fixture
def spades_1_maxk_17() -> Tuple[str, int]:
    return os.path.join('tests', 'data', 'test_contigs_spades_1.fasta.gz'), 17
# end def spades_1_maxk_17


# === Fixtures for function `src.contigs.assign_multiplty` ===

@pytest.fixture
def contig_collection_spades_0() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
        25 # maxk
    )

    return contig_collection
# end def contig_collection_spades_0

@pytest.fixture
def contig_collection_spades_1() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_1.fasta.gz'),
        17 # maxk
    )

    return contig_collection
# end def contig_collection_spades_1

@pytest.fixture
def contig_collection_a5_0() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta'),
        25 # maxk
    )

    return contig_collection
# end def contig_collection_a5_0

@pytest.fixture
def contig_collection_mix_0() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_mix_0.fasta'),
        25 # maxk
    )

    return contig_collection
# end def contig_collection_mix_0

@pytest.fixture
def contig_collection_spades_no_multplty_0() -> cnt.ContigCollection:
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data', 'test_contigs_spades_no_multplty_0.fasta'),
        25 # maxk
    )

    return contig_collection
# end def contig_collection_spades_no_multplty_0



# === Test classes ===

class TestRC:
    # Class for testing function `src.contigs._rc`

    def test_rc_some_sequence(self, some_sequence):
        # Test reverse-complement of a simple sequence
        expected: str = 'AGCCAGATTCATTATGGT'

        assert cnt._rc(some_sequence) == expected
    # end def test_rc_some_sequence

    def test_rc_degenerate_sequence(self, degenerate_sequence):
        # Test reverse-complement of a degenerate sequence
        expected: str = 'ANCCHGATTCATTSTYGT'

        assert cnt._rc(degenerate_sequence) == expected
    # end def test_rc_degenerate_sequence
# end class TestRC


class TestValidateFasta:
    # Class for testing function `src.contigs._validate_fasta`

    def test_validate_fasta_ok(self, valid_fasta: Tuple[str, str]):
        # Test validateion of a valid fasta record
        try:
            cnt._validate_fasta(valid_fasta[0], valid_fasta[1])
        except ValueError:
            pytest.fail('Test `test_validate_fasta_ok` failed!')
        # end try
    # end def test_validate_fasta_ok

    def test_validate_fasta_invalid_header(self, invalid_fasta_header: Tuple[str, str]):
        # Test validateion of a fasta record with invalid header
        with pytest.raises(ValueError):
            cnt._validate_fasta(invalid_fasta_header[0], invalid_fasta_header[1])
        # end with
    # end def test_validate_fasta_invalid_header

    def test_validate_fasta_empty_seq(self, invalid_fasta_seq: Tuple[str, str]):
        # Test validateion of a fasta record with invalid sequence
        with pytest.raises(ValueError):
            cnt._validate_fasta(invalid_fasta_seq[0], invalid_fasta_seq[1])
        # end with
    # end def test_validate_fasta_empty_seq

    def test_validate_fasta_invalid_seq(self, invalid_empty_seq: Tuple[str, str]):
        # Test validateion of a fasta record with empty sequence
        with pytest.raises(ValueError):
            cnt._validate_fasta(invalid_empty_seq[0], invalid_empty_seq[1])
        # end with
    # end def test_validate_fasta_invalid_seq
# end class TestValidateFasta


class TestFastaGeerator:
    # Class for testing generator `src.contigs._fasta_generator`
    def test_fasta_generator_plain(self, plain_text_fasta: str):
        # Check how the generator parses a plain fasta file
        expected_collection: Sequence[Tuple[str, str]] = tuple([
            (
                'NODE_1_length_938688_cov_24.261358',
                'ACCATAATGAATCTGGCTATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCGAATTTTGTGGCACT'
            ),
            (
                'NODE_2_length_622710_cov_18.276868',
                'CGAATTTTGTGGCACTGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCT'
            ),
            (
                'NODE_3_length_396034_cov_43.731463',
                'TTTCACAGTATTCACTGAGACGCTATAACAATACTAGATGGAATTTCACAGTATTCACTGAGAC'
            ),
            (
                'NODE_4_length_381490_cov_28.678435',
                'AATACCTACAACTTGTGCTAATGACCCTGTGGGTTTTACACTTAAAAACACAG'
            ),
        ])

        expected: Tuple[str, str]
        obtained: Tuple[str, str]
        for expected, obtained in zip(expected_collection, cnt._fasta_generator(plain_text_fasta)):
            assert obtained == expected
        # end for
    # end def test_fasta_generator_plain

    def test_fasta_generator_gzipped(self, gzipped_fasta: str):
        # Check how the generator parses a gzipped fasta file
        expected_collection: Sequence[Tuple[str, str]] = tuple([
            (
            'NODE_1_length_70_cov_24.261358',
            'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA',
            ),
            (
            'NODE_2_length_70_cov_24.261358',
            'TCTGTTCTCTAAACTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATTAAAGGTTTA',
            ),
            (
            'NODE_3_length_70_cov_24.261358',
            'TAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTG',
            ),
            (
            'NODE_4_length_70_cov_24.261358',
            'TTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTACACGGACGAAACCGTAA',
            ),
            (
            'NODE_5_length_70_cov_24.261358',
            'CCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTAC',
            ),
            (
            'NODE_6_length_70_cov_24.261358',
            'TCGTTGAAACCAGGTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGG',
            ),
            (
            'NODE_7_length_50_cov_24.261358',
            'CTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACACCTGGTTT',
            ),
        ])

        expected: Tuple[str, str]
        obtained: Tuple[str, str]
        for expected, obtained in zip(expected_collection, cnt._fasta_generator(gzipped_fasta)):
            assert obtained == expected
        # end for
    # end def test_fasta_generator_gzipped
# end class TestFastaGeerator


class TestIsSpadesHeader:
    # Class for testing function `src.contigs._is_spades_header`
    def test_is_spades_header_valid_spades(self, valid_spades_header: str):
        # Test detection of valid SPAdes header
        assert cnt._is_spades_header(valid_spades_header) == True
    # end def test_is_spades_header_valid_spades

    def test_is_spades_header_cryptic_invalid_spades(self, cryptic_invalid_spades_header: str):
        # Test detection of cryptic invalid SPAdes header.
        # It is invalid, but coverage can't be parsed.
        assert cnt._is_spades_header(cryptic_invalid_spades_header) == True
    # end def test_is_spades_header_cryptic_invalid_spades

    def test_is_spades_header_valid_a5(self, valid_a5_header: str):
        # Test detection of valid a5 header
        assert cnt._is_spades_header(valid_a5_header) == False
    # end def test_is_spades_header_valid_a5
# end class TestIsSpadesHeader


class TestFormatContigName:
    # Class for testing function `src.contigs._format_contig_name`
    def test_format_contig_name_valid_spades(self, valid_spades_header: str):
        # Test simplification of valid SPAdes header
        expected: str = 'NODE_1'
        assert cnt._format_contig_name(valid_spades_header) == expected
    # end def test_is_spades_header_valid_spades

    def test_format_contig_name_cryptic_invalid_spades(self, cryptic_invalid_spades_header: str):
        # Test simplification of cryptic invalid SPAdes header.
        # It is invalid, but coverage can't be parsed.
        expected: str = 'NODE_1'
        assert cnt._format_contig_name(cryptic_invalid_spades_header) == expected
    # end def test_is_spades_header_cryptic_invalid_spades

    def test_format_contig_name_valid_a5(self, valid_a5_header: str):
        # Test simplification of valid a5 header
        expected: str = 'scaffold_24'
        assert cnt._format_contig_name(valid_a5_header) == expected
    # end def test_format_contig_name_valid_a5
# end class TestFormatContigName


class TestParseCoverage:
    # Class for testing function `src.contigs._format_contig_name`
    def test_parse_coverage_valid_spades(self, valid_spades_header: str):
        # Test simplification of valid SPAdes header
        expected: float = 24.26
        assert cnt._parse_coverage(valid_spades_header) == expected
    # end def test_parse_coverage_valid_spades

    def test_parse_coverage_cryptic_invalid_spades(self, cryptic_invalid_spades_header: str):
        # Test simplification of cryptic invalid SPAdes header.
        # It is invalid, but coverage can't be parsed.
        expected: float = None
        assert cnt._parse_coverage(cryptic_invalid_spades_header) == expected
    # end def test_parse_coverage_cryptic_invalid_spades

    def test_parse_coverage_valid_a5(self, valid_a5_header: str):
        # Test simplification of valid a5 header
        expected: float = None
        assert cnt._parse_coverage(valid_a5_header) == expected
    # end def test_parse_coverage_valid_a5
# end class TestParseCoverage


class TestCalcGcContent:
    # Class for testing function `src.contigs._calc_gc_content`
    def test_calc_gc_content_some_seq(self, some_sequence):
        expected: float = 38.89
        assert abs(cnt._calc_gc_content(some_sequence) - expected) < 1e-2
    # end def test_calc_gc_content_some_seq

    def test_calc_gc_content_degenerate_seq(self, degenerate_sequence):
        expected: float = 33.33
        assert abs(cnt._calc_gc_content(degenerate_sequence) - expected) < 1e-2
    # end def test_calc_gc_content_degenerate_seq
# end class TestCalcGcContent


class TestGetContigCollection:
    # Class for testing function `src.contigs.get_contig_collection`

    GetAttrFunc = Callable[[cnt.Contig], Any]

    @staticmethod
    def _attribute_generator() -> Generator[Tuple[str, GetAttrFunc, type], None, None]:
        # Generator yields stuff for testing attributes of Class `src.contigs.Contig`:
        #  1. Name of and attribute.
        #  2. Function for accessing this attribute.
        #  3. Type to which an attribure must belong.

        yield 'name', lambda x: x.name, str
        yield 'length', lambda x: x.length, int
        yield 'cov', lambda x: x.cov, float
        yield 'gc_content', lambda x: x.gc_content, float
        yield 'start', lambda x: x.start, str
        yield 'rcstart', lambda x: x.rcstart, str
        yield 'end', lambda x: x.end, str
        yield 'rcend', lambda x: x.rcend, str
    # end def _test_attribute

    def test_get_contig_collection_spades_1(self, spades_1_maxk_17):
        # Function for testing function `src.contigs.get_contig_collection`

        infpath: str = spades_1_maxk_17[0]
        maxk: int = spades_1_maxk_17[1]

        contig_collection: cnt.ContigCollection = cnt.get_contig_collection(infpath, maxk)

        # Check length of the collection
        expected_len: int = 7
        assert len(contig_collection) == expected_len

        # Check if all items are of class `src.contigs.Contig`
        assert all(map(lambda x: isinstance(x, cnt.Contig), contig_collection))

        # Check if all `Contig`s' attributes are not None (except `multplty`) and of proper types
        attr_name: str
        get_attr: GetAttrFunc
        attr_type: type
        for attr_name, get_attr, attr_type in self._attribute_generator():
            try:
                # Check if it is not None
                assert all(map(lambda x: not get_attr(x) is None, contig_collection))
                # Check type
                assert all(map(lambda x: isinstance(get_attr(x), attr_type), contig_collection))
            except AttributeError:
                pytest.fail('Error: Contig instance attribute has no `{}` attribute.'\
                    .format(attr_name))
            # end try
        # end for

    # end def test_get_contig_collection_spades_1
# end class TestGetContigCollection


class TestAssignMultiplty:
    # Class for testing function `src.contigs.assign_multiplty`

    def test_assign_multiplty_spades_0(self, contig_collection_spades_0):
        # Test how the function assigns multiplicity to contigs from
        #   file `test_contigs_spades_0.fasta

        # "Rename" variable
        contig_collection: cnt.ContigCollection = contig_collection_spades_0

        # Assign
        cnt.assign_multiplty(contig_collection)

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
        contig_collection: cnt.ContigCollection = contig_collection_spades_1

        # Assign
        cnt.assign_multiplty(contig_collection)

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
        contig_collection: cnt.ContigCollection = contig_collection_a5_0

        # Assign
        cnt.assign_multiplty(contig_collection)

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
        contig_collection: cnt.ContigCollection = contig_collection_mix_0

        # Assign
        cnt.assign_multiplty(contig_collection)

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
        contig_collection: cnt.ContigCollection = contig_collection_spades_no_multplty_0

        # Assign
        cnt.assign_multiplty(contig_collection)

        expected_multplts: Sequence[float] = tuple([1.0] * 4)

        obtained_multplts: Sequence[float] = tuple(map(lambda x: x.multplty, contig_collection))

        # Compare lengths
        assert len(obtained_multplts) == len(expected_multplts)

        # Compare multiplicities
        for expected, obtained in zip(expected_multplts, obtained_multplts):
            assert abs(obtained - expected) < 1e-1
        # end for
    # end def test_assign_multiplty_spades_no_multplty_0
# end class TestAssignMultiplty
