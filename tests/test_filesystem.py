# -*- coding: utf-8 -*-

import os
import pytest

import src.filesystem as fls

# === Fixtures for `src.filesystem.is_fasta` ===
# (some also for `src.filesystem._bname_no_fasta_ext` and `src.filesystem._make_numbered_prefix`)

@pytest.fixture(scope='session')
def fpath_fasta():
    # Returns path to file having `.fasta` extention
    return os.path.join('some', 'happy', 'file.fasta')
# end def fpath_fasta


@pytest.fixture
def fpath_fa():
    # Returns path to file having `.fa` extention
    return os.path.join('some', 'happy', 'file.fa')
# end def fpath_fa


@pytest.fixture
def fpath_fasta_gz():
    # Returns path to file having `.fasta.gz` extention
    return os.path.join('some', 'happy', 'file.fasta.gz')
# end def fpath_fasta_gz


@pytest.fixture
def fpath_fa_gz():
    # Returns path to file having `.fa.gz` extention
    return os.path.join('some', 'happy', 'file.fa.gz')
# end def fpath_fa_gz


@pytest.fixture
def fpath_fasta_bz2():
    # Returns path to file having `.fasta.bz2` extention
    return os.path.join('some', 'happy', 'file.fasta.bz2')
# end def fpath_fasta_bz2


@pytest.fixture
def fpath_fa_bz2():
    # Returns path to file having `.fa.bz2` extention
    return os.path.join('some', 'happy', 'file.fa.bz2')
# end def fpath_fa_bz2


@pytest.fixture
def fpath_fast():
    # Returns path to file having `.fast` extention
    return os.path.join('some', 'happy', 'file.fast')
# end def fpath_fast


@pytest.fixture
def fpath_txt():
    # Returns path to file having `.txt` extention
    return os.path.join('some', 'happy', 'file.txt')
# end def fpath_txt


# === Fixtures for `src.filesystem.conf_prefix` ===

@pytest.fixture(scope='session')
def empty_outdir(tmpdir_factory):
    # Empty output directory
    outdpath: str = tmpdir_factory.mktemp('empty-outdir')
    fls.make_outdir(outdpath)
    return outdpath
# end def empty_outdir


@pytest.fixture(scope='session')
def outdir_single_file(tmpdir_factory, fpath_fasta):
    # Output directory with single `file__combinator_adjacent_contigs.tsv` file
    outdpath: str = tmpdir_factory.mktemp('outdir-single-file')
    fls.make_outdir(outdpath)

    # Configure path to origin output file
    origin_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))
    )

    # Create origin output file
    with open(origin_fpath, 'w') as tmpfile:
        pass
    # end with
    return outdpath
# end def outdir_single_file


@pytest.fixture(scope='session')
def outdir_two_files(tmpdir_factory, fpath_fasta):
    # Output directory with following files:
    #   `file_combinator_adjacent_contigs.tsv`
    #   `file.1_combinator_adjacent_contigs.tsv`
    outdpath: str = tmpdir_factory.mktemp('outdir-two-files')
    fls.make_outdir(outdpath)

    # Configure path to origin output file
    origin_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))
    )

    # Create origin output file
    with open(origin_fpath, 'w') as tmpfile:
        pass
    # end with

    # Create numbered (number 1) file
    with open(
        os.path.join(
            outdpath,
            '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))),
        'w') as tmpfile:
        pass
    # end with
    return outdpath
# end def outdir_two_files


@pytest.fixture(scope='session')
def outdir_single_numbered(tmpdir_factory, fpath_fasta):
    # Output directory with following files:
    #   `file.1_combinator_adjacent_contigs.tsv`
    outdpath: str = tmpdir_factory.mktemp('outdir-two-files')
    fls.make_outdir(outdpath)

    # Configure path to origin output file
    origin_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))
    )

    # Create origin output file
    with open(origin_fpath, 'w') as tmpfile:
        pass
    # end with

    # Create numbered (number 1) file
    with open(
        os.path.join(
            outdpath,
            '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))),
        'w') as tmpfile:
        pass
    # end with

    # And remove origin file
    os.unlink(origin_fpath)

    return outdpath
# end def outdir_single_numbered


@pytest.fixture(scope='session')
def outdir_three_files(tmpdir_factory, fpath_fasta):
    # Output directory with following files:
    #   `file_combinator_adjacent_contigs.tsv`
    #   `file.1_combinator_adjacent_contigs.tsv`
    #   `file.2_combinator_adjacent_contigs.tsv`
    outdpath: str = tmpdir_factory.mktemp('outdir-three-files')
    fls.make_outdir(outdpath)

    # Configure path to origin output file
    origin_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))
    )

    # Create origin output file
    with open(origin_fpath, 'w') as tmpfile:
        pass
    # end with

    # Create numbered (numbers 1, 2) files
    i: int
    for i in range(2):
        with open(
            os.path.join(
                outdpath,
                '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))),
            'w') as tmpfile:
            pass
        # end with
    # end for
    return outdpath
# end def outdir_three_files


@pytest.fixture(scope='session')
def outdir_two_files_gap(tmpdir_factory, fpath_fasta):
    # Output directory with following files:
    #   `file_combinator_adjacent_contigs.tsv`
    #   `file.2_combinator_adjacent_contigs.tsv`
    outdpath: str = tmpdir_factory.mktemp('outdir-three-files')
    fls.make_outdir(outdpath)

    # Configure path to origin output file
    origin_fpath: str = os.path.join(
        outdpath,
        '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))
    )

    # Create origin output file
    with open(origin_fpath, 'w') as tmpfile:
        pass
    # end with

    # Create numbered (numbers 1, 2) files
    i: int
    for i in range(2):
        with open(
            os.path.join(
                outdpath,
                '{}_combinator_adjacent_contigs.tsv'.format(fls.conf_prefix(fpath_fasta, outdpath))),
            'w') as tmpfile:
            pass
        # end with
    # end for

    # And remove file number 1
    os.unlink(os.path.join(outdpath, 'file.1_combinator_adjacent_contigs.tsv'))

    return outdpath
# end def outdir_two_files_gap


# === Test classes ===

class TestIsFasta:
    # Class for testing `src.filesystem.is_fasta`

    def test_is_fasta_fasta(self, fpath_fasta):
        # Should recognize file with `.fasta` extention
        assert fls.is_fasta(fpath_fasta) == True
    # end def test_is_fasta_fasta

    def test_is_fasta_fa(self, fpath_fa):
        # Should recognize file with `.fa` extention
        assert fls.is_fasta(fpath_fa) == True
    # end def test_is_fasta_fa

    def test_is_fasta_fasta_gz(self, fpath_fasta_gz):
        # Should recognize file with `.fasta.gz` extention
        assert fls.is_fasta(fpath_fasta_gz) == True
    # end def test_is_fasta_fasta_gz

    def test_is_fasta_fa_gz(self, fpath_fa_gz):
        # Should recognize file with `.fa.gz` extention
        assert fls.is_fasta(fpath_fa_gz) == True
    # end def test_is_fasta_fa_gz

    def test_is_fasta_fasta_bz2(self, fpath_fasta_bz2):
        # Should not recognize file with `.fasta.bz2` extention
        assert fls.is_fasta(fpath_fasta_bz2) == False
    # end def test_is_fasta_fasta_bz2

    def test_is_fasta_fa_bz2(self, fpath_fa_bz2):
        # Should not recognize file with `.fa.bz2` extention
        assert fls.is_fasta(fpath_fa_bz2) == False
    # end def test_is_fasta_fa_bz2

    def test_is_fasta_fast(self, fpath_fast):
        # Should not recognize file with `.fast` extention
        assert fls.is_fasta(fpath_fast) == False
    # end def test_is_fasta_fasta

    def test_is_fasta_txt(self, fpath_txt):
        # Should not recognize file with `.txt` extention
        assert fls.is_fasta(fpath_txt) == False
    # end def test_is_fasta_fasta
# end class TestIsFasta


class TestMakeOutdir:
    # Class for testing `src.filesystem.make_outdir`

    def test_make_outdir(self, tmpdir):
        # Should successfully make a directory in `tempdir`
        #   and create a file in it.

        outdpath: str = os.path.join(tmpdir, 'outdir')
        fls.make_outdir(outdpath)

        assert os.path.isdir(outdpath)

        tesf_fpath = os.path.join(outdpath, 'test.file.txt')
        with open(tesf_fpath, 'w') as testfile:
            pass
        # end with

        assert os.path.isfile(tesf_fpath)
    # end def test_make_outdir
# end class TestMakeOutdir


class TestBnameNoFastaExt:
    # Class for testing `src.filesystem._bname_no_fasta_ext`

    def test_bname_no_fasta_ext_fasta(self, fpath_fasta):
        # Should successfully return basename without `.fasta` extention
        assert fls._bname_no_fasta_ext(fpath_fasta) == 'file'
    # end def test_bname_no_fasta_ext_fasta

    def test_bname_no_fasta_ext_fa(self, fpath_fa):
        # Should successfully return basename without `.fa` extention
        assert fls._bname_no_fasta_ext(fpath_fa) == 'file'
    # end def test_bname_no_fasta_ext_fa

    def test_bname_no_fasta_ext_fasta_gz(self, fpath_fasta_gz):
        # Should successfully return basename without `.fasta.gz` extention
        assert fls._bname_no_fasta_ext(fpath_fasta_gz) == 'file'
    # end def test_bname_no_fasta_ext_fasta_gz

    def test_bname_no_fasta_ext_fa_gz(self, fpath_fa_gz):
        # Should successfully return basename without `.fa.gz` extention
        assert fls._bname_no_fasta_ext(fpath_fa_gz) == 'file'
    # end def test_bname_no_fasta_ext_fa_gz
# end class TestBnameNoFastaExt


class TestMakeNumberedPrefix:
    # Class for testing `src.filesystem._make_numbered_prefix`

    def test_make_numbered_prefix_fasta(self, fpath_fasta):
        # Should successfully return prefix (number = 1) without `.fasta` extention
        assert fls._make_numbered_prefix(fpath_fasta, 1) == 'file.1'
    # end def test_make_numbered_prefix_fasta

    def test_make_numbered_prefix_fa(self, fpath_fa):
        # Should successfully return prefix (number = 1) without `.fa` extention
        assert fls._make_numbered_prefix(fpath_fa, 1) == 'file.1'
    # end def test_make_numbered_prefix_fa

    def test_make_numbered_prefix_fasta_gz(self, fpath_fasta_gz):
        # Should successfully return prefix (number = 1) without `.fasta.gz` extention
        assert fls._make_numbered_prefix(fpath_fasta_gz, 1) == 'file.1'
    # end def test_make_numbered_prefix_fasta_gz

    def test_make_numbered_prefix_fa_gz(self, fpath_fa_gz):
        # Should successfully return prefix (number = 1) without `.fa.gz` extention
        assert fls._make_numbered_prefix(fpath_fa_gz, 1) == 'file.1'
    # end def test_make_numbered_prefix_fa_gz
# end class TestMakeNumberedPrefix


class TestConfPrefix:
    # Class for testing `src.filesystem.conf_prefix`

    def test_conf_prefix_empty(self, fpath_fasta, empty_outdir):
        # Prefix should be just `file`
        assert fls.conf_prefix(fpath_fasta, empty_outdir) == 'file'
    # end def test_conf_prefix_empty

    def test_conf_prefix_single(self, fpath_fasta, outdir_single_file):
        # Prefix should be `file.1`
        assert fls.conf_prefix(fpath_fasta, outdir_single_file) == 'file.1'
    # end def test_conf_prefix_single

    def test_conf_prefix_two(self, fpath_fasta, outdir_two_files):
        # Prefix should be `file.2`
        assert fls.conf_prefix(fpath_fasta, outdir_two_files) == 'file.2'
    # end def test_conf_prefix_two

    def test_conf_prefix_single_numbered(self, fpath_fasta, outdir_single_numbered):
        # Prefix should be just `file`
        assert fls.conf_prefix(fpath_fasta, outdir_single_numbered) == 'file'
    # end def test_conf_prefix_single_numbered

    def test_conf_prefix_three(self, fpath_fasta, outdir_three_files):
        # Prefix should be `file.3`
        assert fls.conf_prefix(fpath_fasta, outdir_three_files) == 'file.3'
    # end def test_conf_prefix_three
# end class TestConfPrefix
