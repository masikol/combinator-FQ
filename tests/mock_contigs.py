# -*- coding: utf-8 -*-

import os
import pytest
from typing import Tuple

import src.contigs as cnt
import src.overlaps as ovl




MockContigsFixture = Tuple[cnt.ContigCollection, ovl.OverlapCollection]


@pytest.fixture(scope='session')
def mock_contigs_spades_0() -> MockContigsFixture:

    contig_collection, overlap_collection = _do_combinator_work(
        infpath=os.path.join('tests', 'data/test_contigs_spades_0.fasta'),
        mink=16,
        maxk=25
    )
    return contig_collection, overlap_collection
# end def mock_contigs_spades_0


@pytest.fixture(scope='session')
def mock_contigs_spades_1() -> MockContigsFixture:

    contig_collection, overlap_collection = _do_combinator_work(
        infpath=os.path.join('tests', 'data/test_contigs_spades_1.fasta'),
        mink=8,
        maxk=17
    )
    return contig_collection, overlap_collection
# end def mock_contigs_spades_1


@pytest.fixture(scope='session')
def mock_contigs_a5_0() -> MockContigsFixture:

    contig_collection, overlap_collection = _do_combinator_work(
        infpath=os.path.join('tests', 'data/test_contigs_a5_0.fasta'),
        mink=16,
        maxk=25
    )
    return contig_collection, overlap_collection
# end def mock_contigs_a5_0


@pytest.fixture(scope='session')
def mock_contigs_mix_0() -> MockContigsFixture:

    contig_collection, overlap_collection = _do_combinator_work(
        infpath=os.path.join('tests', 'data/test_contigs_mix_0.fasta'),
        mink=16,
        maxk=25
    )
    return contig_collection, overlap_collection
# end def mock_contigs_mix_0


@pytest.fixture(scope='session')
def mock_contigs_spades_no_multplty_0() -> MockContigsFixture:

    contig_collection, overlap_collection = _do_combinator_work(
        infpath=os.path.join('tests', 'data/test_contigs_spades_no_multplty_0.fasta'),
        mink=16,
        maxk=25
    )
    return contig_collection, overlap_collection
# end def mock_contigs_spades_no_multplty_0


def _do_combinator_work(infpath: str, mink: int, maxk: int) -> MockContigsFixture:

    # Read contigs
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        infpath,
        maxk
    )
    # Assign multiplicity to contigs
    cnt.assign_multiplty(contig_collection)

    # Detect adjacent contigs
    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def _do_combinator_work
