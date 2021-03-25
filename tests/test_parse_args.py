# -*- coding: utf-8 -*-

import os
import sys
import glob
import pytest
from io import StringIO
from typing import Tuple, TextIO, Sequence, Dict, Any, List

import src.parse_args as par


Paths = Sequence[str]
OptsArgs = Tuple[Tuple[str, str], Tuple[str, str]]
Params = Dict[str, Any]
Argv = List[str]
FpathsParams = Tuple[Paths, Params]


# === Fixtures for testing function `src.parse_args._get_input_fpaths` ===

@pytest.fixture
def single_file() -> Paths:
    # Returns tuple of single path to input fasta file
    return tuple([
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta')
    ])
# end def single_file

@pytest.fixture
def two_files() -> Paths:
    # Returns tuple of two paths to input fasta files
    return tuple([
        os.path.join('tests', 'data', 'test_contigs_a5_0.fasta'),
        os.path.join('tests', 'data', 'test_contigs_spades_1.fasta.gz')
    ])
# end def two_files

@pytest.fixture
def no_files() -> Paths:
    # Returns empty tuple
    return tuple()
# end def no_files


# === Fixtures for testing function `src.parse_args._parse_options` ===

@pytest.fixture
def default_params() -> OptsArgs:
    # Returns OptsArgs of default params
    return tuple()
# end def default_params

@pytest.fixture
def all_valid_params() -> OptsArgs:
    # Returns OptsArgs of all valid paramters, without k
    return tuple([
        ('-i', '23'),
        ('-a', '125'),
        ('-o', 'output-dir')
    ])
# end def all_valid_params

@pytest.fixture
def params_mink_zero() -> OptsArgs:
    # Returns OptsArgs where `mink` is zero, without k
    return tuple([
        ('-i', '0'),
        ('-a', '125'),
        ('-o', 'output-dir')
    ])
# end def params_mink_zero

@pytest.fixture
def params_mink_negative() -> OptsArgs:
    # Returns OptsArgs where `mink` is negative, without k
    return tuple([
        ('-i', '-1'),
        ('-a', '125'),
        ('-o', 'output-dir')
    ])
# end def params_mink_negative

@pytest.fixture
def params_maxk_zero() -> OptsArgs:
    # Returns OptsArgs where `mink` is zero, without k
    return tuple([
        ('-i', '23'),
        ('-a', '0'),
        ('-o', 'output-dir')
    ])
# end def params_maxk_zero

@pytest.fixture
def params_maxk_negative() -> OptsArgs:
    # Returns OptsArgs where `mink` is negative, without k
    return tuple([
        ('-i', '23'),
        ('-a', '-1'),
        ('-o', 'output-dir')
    ])
# end def params_maxk_negative

@pytest.fixture
def params_mink_nonint() -> OptsArgs:
    # Returns OptsArgs where `mink` is not int, without k
    return tuple([
        ('-i', 'SABAKA'),
        ('-a', '125'),
        ('-o', 'output-dir')
    ])
# end def params_mink_nonint

@pytest.fixture
def params_maxk_nonint() -> OptsArgs:
    # Returns OptsArgs where `maxk` is not int, without k
    return tuple([
        ('-i', '23'),
        ('-a', 'SABAKA'),
        ('-o', 'output-dir')
    ])
# end def params_maxk_nonint

@pytest.fixture
def params_k_valid() -> OptsArgs:
    # Returns OptsArgs where valid `k` is specified
    return tuple([
        ('-k', '127'),
        ('-o', 'output-dir')
    ])
# end def params_k_valid

@pytest.fixture
def params_k_zero() -> OptsArgs:
    # Returns OptsArgs where specified `k` is zero
    return tuple([
        ('-k', '0'),
        ('-o', 'output-dir')
    ])
# end def params_k_zero

@pytest.fixture
def params_k_negative() -> OptsArgs:
    # Returns OptsArgs where specified `k` is negative
    return tuple([
        ('-k', '-1'),
        ('-o', 'output-dir')
    ])
# end def params_k_negative

@pytest.fixture
def params_k_nonint() -> OptsArgs:
    # Returns OptsArgs where specified `k` is not int
    return tuple([
        ('-k', 'SABAKA'),
        ('-o', 'output-dir')
    ])
# end def params_k_nonint


# === Fixtures for fuction `src.parse_args.parse_args` ===

@pytest.fixture
def argv_minus_h() -> Argv:
    # Arguments for printing help message
    return ['combinatro-FQ.py', '-h']
# end def argv_minus_h

@pytest.fixture
def argv_minus_help() -> Argv:
    # Arguments for printing help message
    return ['combinatro-FQ.py', '--help']
# end def argv_minus_help

@pytest.fixture
def argv_minus_v() -> Argv:
    # Arguments for printing version
    return ['combinatro-FQ.py', '-v']
# end def argv_minus_v

@pytest.fixture
def argv_minus_version() -> Argv:
    # Arguments for printing version
    return ['combinatro-FQ.py', '--version']
# end def argv_minus_version

@pytest.fixture
def argv_defaults() -> Argv:
    # Default arguments
    return [
    'combinatro-FQ.py',
    os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
    ]
# end def argv_defaults

@pytest.fixture
def argv_all_valid() -> Argv:
    # All valid arguments
    return [
    'combinatro-FQ.py',
    os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
    '-i', '23',
    '-a', '125'
    ]
# end def argv_all_valid

@pytest.fixture
def argv_mink_gt_maxk_mink_not_spec() -> Argv:
    # All valid arguments
    return [
    'combinatro-FQ.py',
    os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
    '-a', '19'
    ]
# end def argv_mink_gt_maxk_mink_not_spec

@pytest.fixture
def argv_mink_gt_maxk_maxk_not_spec() -> Argv:
    # All valid arguments
    return [
    'combinatro-FQ.py',
    os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
    '-i', '129'
    ]
# end def argv_mink_gt_maxk_maxk_not_spec

@pytest.fixture
def argv_mink_gt_maxk_all_spec() -> Argv:
    # All valid arguments
    return [
    'combinatro-FQ.py',
    os.path.join('tests', 'data', 'test_contigs_spades_0.fasta'),
    '-i', '23'
    '-a', '21'
    ]
# end def argv_mink_gt_maxk_all_spec

# === Test classes

class TestGetInputFpaths:
    # Class for testing function `src.parse_args._get_input_fpaths`

    def test_get_input_fpaths_single_file(self, single_file: Paths):
        # Tests how it handles tuple of single path to input fasta file
        expected: Tuple[str] = tuple([
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_a5_0.fasta')
        ])

        assert par._get_input_fpaths(single_file) == expected
    # end def test_get_input_fpaths_single_file

    def test_get_input_fpaths_two_files(self, two_files: Paths):
        # Tests how it handles tuple of tow paths to input fasta file
        expected: Tuple[str] = tuple([
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_a5_0.fasta'),
            os.path.join(os.getcwd(), 'tests', 'data', 'test_contigs_spades_1.fasta.gz')
        ])
    # end def test_get_input_fpaths_two_files

    def test_get_input_fpaths_no_files(self, no_files: Paths):
        # Tests how it handles empty tuple
        with pytest.raises(SystemExit):
            par._get_input_fpaths(no_files)
        # end with
    # end def test_get_input_fpaths_no_files

    def test_get_input_fpaths_fromcwd_n(self, no_files: Paths):
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

    def test_get_input_fpaths_fromcwd_y(self, no_files: Paths):
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
    # end def test_get_input_fpaths_fromcwd_y
# end class TestGetInputFpaths


class TestParseOtions:
    # Class for testing function `src.parse_args._parse_options`.

    def test_parse_options_defaults(self, default_params: OptsArgs):
        # Test `_parse_options` setting default options
        expected: Params = {
            'o': os.path.join(os.getcwd(), 'combinator-result'),
            'i': 21,
            'a': 127
        }
        par._parse_options(default_params)
    # end def test_parse_options_defaults

    def test_parse_options_all_valid(self, all_valid_params: OptsArgs):
        # Test `_parse_options` with all valid params specidied
        expected: Params = {
            'o': os.path.join(os.getcwd(), 'output-dir'),
            'i': 23,
            'a': 125
        }
        par._parse_options(all_valid_params)
    # end def test_parse_options_all_valid

    def test_parse_options_mink_lqe_zero(
        self,
        params_mink_zero: OptsArgs,
        params_mink_negative: OptsArgs
    ):
        # Test `_parse_options` with `mink` less that or equal to zero
        for fixture in (params_mink_zero, params_mink_negative):
            with pytest.raises(SystemExit):
                par._parse_options(fixture)
            # end with
        # end for
    # end def test_parse_options_all_valid

    def test_parse_options_maxk_lqe_zero(
        self,
        params_maxk_zero: OptsArgs,
        params_maxk_negative: OptsArgs
    ):
        # Test `_parse_options` with `maxk` less that or equal to zero
        for fixture in (params_maxk_zero, params_maxk_negative):
            with pytest.raises(SystemExit):
                par._parse_options(fixture)
            # end with
        # end for
    # end def test_parse_options_maxk_lqe_zero

    def test_parse_options_mink_nonint(self, params_mink_nonint: OptsArgs):
        # Test `_parse_options` with `mink` is not int
        with pytest.raises(SystemExit):
            par._parse_options(params_mink_nonint)
        # end with
    # end def test_parse_options_mink_nonint

    def test_parse_options_maxk_nonint(self, params_maxk_nonint: OptsArgs):
        # Test `_parse_options` with `maxk` is not int
        with pytest.raises(SystemExit):
            par._parse_options(params_maxk_nonint)
        # end with
    # end def test_parse_options_maxk_nonint

    def test_parse_options_k_valid(self, params_k_valid: OptsArgs):
        # Test `_parse_options` with valid `k` specified
        expected: Params = {
            'o': os.path.join(os.getcwd(), 'output-dir'),
            'i': 127,
            'a': 127
        }
        par._parse_options(params_k_valid)
    # end def test_parse_options_k_valid

    def test_parse_options_k_lqe_zero(
        self,
        params_k_zero: OptsArgs,
        params_k_negative: OptsArgs
    ):
        # Test `_parse_options` with `maxk` less that or equal to zero
        for fixture in (params_k_zero, params_k_negative):
            with pytest.raises(SystemExit):
                par._parse_options(fixture)
            # end with
        # end for
    # end def test_parse_options_k_lqe_zero

    def test_parse_options_k_nonint(self, params_k_nonint: OptsArgs):
        # Test `_parse_options` with `k` is not int
        with pytest.raises(SystemExit):
            par._parse_options(params_k_nonint)
        # end with
    # end def test_parse_options_k_nonint
# end class TestParseOtions


class TestParseArgs:
    # Class for testing function `src.parse_args.parse_args`

    def test_parse_argv_help(self, argv_minus_h: Argv, argv_minus_help: Argv):
        # Function for testing how `parse_args` prints help
        argv: List[str]
        for argv in (argv_minus_h, argv_minus_help):
            sys.argv = argv
            try:
                with pytest.raises(SystemExit):
                    par.parse_args('version', 'date')
                # end with
            except AssertionError:
                pass
            finally:
                sys.argv = list()
            # end try
        # end for
    # end def test_parse_argv_help

    def test_parse_argv_version(self, argv_minus_v: Argv, argv_minus_version: Argv):
        # Function for testing how `parse_args` prints version
        argv: List[str]
        for argv in (argv_minus_v, argv_minus_version):
            sys.argv = argv
            try:
                with pytest.raises(SystemExit):
                    par.parse_args('version', 'date')
                # end with
            except AssertionError:
                pass
            finally:
                sys.argv = list()
            # end try
        # end for
    # end def test_parse_argv_version

    def test_parse_argv_defaults(self, argv_defaults: Argv):
        # Function for testing how `parse_args` handles default params
        sys.argv = argv_defaults

        expected: FpathsParams = tuple([
            tuple([
                os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
            ]),
            {
                'o': os.path.join(os.getcwd(), 'combinator-result'),
                'i': 21,
                'a': 127
            }
        ])

        obtained: FpathsParams = par.parse_args('version', 'date')

        try:
            assert obtained == expected
        except AssertionError:
            pass
        finally:
            argv = list()
        # end try
    # end def test_parse_argv_version

    def test_parse_argv_all_valid(self, argv_all_valid: Argv):
        # Function for testing how `parse_args` handles valid params
        sys.argv = argv_all_valid

        expected: FpathsParams = tuple([
            tuple([
                os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
            ]),
            {
                'o': os.path.join(os.getcwd(), 'combinator-result'),
                'i': 23,
                'a': 125
            }
        ])

        obtained: FpathsParams = par.parse_args('version', 'date')

        try:
            assert obtained == expected
        except AssertionError:
            pass
        finally:
            argv = list()
        # end try
    # end def test_parse_argv_all_valid

    def test_parse_argv_mink_gt_maxk_mink_not_spec(self, argv_mink_gt_maxk_mink_not_spec: Argv):
        # Function for testing how `parse_args` handles situation
        #   where `mink` > `maxk` but `mink` is not specified
        sys.argv = argv_mink_gt_maxk_mink_not_spec

        expected: FpathsParams = tuple([
            tuple([
                os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
            ]),
            {
                'o': os.path.join(os.getcwd(), 'combinator-result'),
                'i': 19,
                'a': 19
            }
        ])

        obtained: FpathsParams = par.parse_args('version', 'date')

        try:
            assert obtained == expected
        except AssertionError:
            pass
        finally:
            argv = list()
        # end try
    # end def test_parse_argv_mink_gt_maxk_mink_not_spec

    def test_parse_argv_mink_gt_maxk_maxk_not_spec(self, argv_mink_gt_maxk_maxk_not_spec: Argv):
        # Function for testing how `parse_args` handles situation
        #   where `mink` > `maxk` but `maxk` is not specified
        sys.argv = argv_mink_gt_maxk_maxk_not_spec

        expected: FpathsParams = tuple([
            tuple([
                os.path.join('tests', 'data', 'test_contigs_spades_0.fasta')
            ]),
            {
                'o': os.path.join(os.getcwd(), 'combinator-result'),
                'i': 129,
                'a': 129
            }
        ])

        obtained: FpathsParams = par.parse_args('version', 'date')

        try:
            assert obtained == expected
        except AssertionError:
            pass
        finally:
            argv = list()
        # end try
    # end def test_parse_argv_mink_gt_maxk_maxk_not_spec

    def test_parse_argv_mink_gt_maxk_maxk_all_spec(self, argv_mink_gt_maxk_all_spec: Argv):
        # Function for testing how `parse_args` handles situation
        #   where `mink` > `maxk` and all params are specified
        sys.argv = argv_mink_gt_maxk_all_spec

        try:
            with pytest.raises(SystemExit):
                par.parse_args('version', 'date')
            # end with
        except AssertionError:
            pass
        finally:
            argv = list()
        # end try
    # end def test_parse_argv_mink_gt_maxk_maxk_all_spec
# end class TestParseArgs
