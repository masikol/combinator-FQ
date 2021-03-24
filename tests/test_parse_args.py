# -*- coding: utf-8 -*-

import os
import sys
import glob
import pytest
from io import StringIO
from typing import Tuple, TextIO

import src.parse_args as par

# === Fixtures ===

@pytest.fixture
def single_file():
    # Returns tuple of single path to input fasta file
    return tuple([
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta')
    ])
# end def single_file

@pytest.fixture
def two_files():
    # Returns tuple of two paths to input fasta files
    return tuple([
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta'),
        os.path.join('tests', 'data', 'test_contigs_spades_1.fasta')
    ])
# end def two_files

@pytest.fixture
def no_files():
    # Returns empty tuple
    return tuple()
# end def no_files

# === Test classes

class TestGetInputFpaths:
    # Class for testing function `src.parse_args._get_input_fpaths`

    def test_get_input_fpaths_single_file(self, single_file):
        # Tests how it handles tuple of single path to input fasta file
        expected: Tuple[str] = tuple([
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_a5_0.fasta')
        ])

        assert par._get_input_fpaths(single_file) == expected
    # end def test_get_input_fpaths_single_file

    def test_get_input_fpaths_two_files(self, two_files):
        # Tests how it handles tuple of tow paths to input fasta file
        expected: Tuple[str] = tuple([
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_a5_0.fasta'),
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_spades_1.fasta')
        ])
    # end def test_get_input_fpaths_two_files

    def test_get_input_fpaths_no_files(self, no_files):
        # Tests how it handles empty tuple
        with pytest.raises(SystemExit):
            par._get_input_fpaths(no_files)
        # end with
    # end def test_get_input_fpaths_no_files

    def test_get_input_fpaths_fromcwd_n(self, no_files):
        # Tests how it handles empty tuple, but there are fasta files in the working dir.
        # However a user enters 'n' or 'N', which terminated the program.
        save_stdin: TextIO = sys.stdin
        save_dir: str = os.getcwd()

        try:
            reply: str
            for reply in ('n', 'N'):
                sys.stdin = StringIO(reply)
                with pytest.raises(SystemExit):
                    par._get_input_fpaths(no_files)
                # end with
            # end for
        except AssertionError:
            pass
        finally:
            sys.stdin = save_stdin
            os.chdir(save_dir)
        # end try
    # end def test_get_input_fpaths_fromcwd_n

    def test_get_input_fpaths_fromcwd_y(self, no_files):
        # Tests how it handles empty tuple, but there are fasta files in the working dir.
        # A user enters 'y' or 'Y', and we get all fasta files in the working dir.
        save_stdin: TextIO = sys.stdin
        save_dir: str = os.getcwd()

        expected: Tuple[str] = glob.glob(os.path.join(os.getcwd(), 'tests', 'data', '*.fasta'))
        os.chdir(os.path.join(os.getcwd(), 'tests', 'data'))

        try:
            reply: str
            for reply in ('y', 'Y'):
                sys.stdin = StringIO(reply)
                assert set(par._get_input_fpaths(no_files)) == set(expected)
            # end for
        except AssertionError:
            pass
        finally:
            sys.stdin = save_stdin
            os.chdir(save_dir)
        # end try

        sys.stdin = save_stdin
        os.chdir(save_dir)
    # end def test_get_input_fpaths_fromcwd_y
# end class TestGetInputFpaths


# class TestParseOtions:
#     # Class for testing function `src.parse_args._parse_options`.

#     def test_parse_options(self, )
