# -*- coding: utf-8 -*-

import os
import pytest

import src.contigs as cnt
import src.overlaps as ovl


@pytest.fixture
def mock_contigs_0():
    mink: int = 16
    maxk: int = 25
    # Read contigs
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data/test_contigs_0.fasta'),
        maxk
    )
    # Assign multiplicity to contigs
    cnt.assign_multiplty(contig_collection)

    # Detect adjacent contigs
    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def mock_contigs_0


@pytest.fixture
def mock_contigs_1():
    mink: int = 8
    maxk: int = 17
    # Read contigs
    contig_collection: cnt.ContigCollection = cnt.get_contig_collection(
        os.path.join('tests', 'data/test_contigs_1.fasta'),
        maxk
    )
    # Assign multiplicity to contigs
    cnt.assign_multiplty(contig_collection)

    # Detect adjacent contigs
    overlap_collection: ovl.OverlapCollection = ovl.detect_adjacent_contigs(
        contig_collection, mink, maxk
    )

    return contig_collection, overlap_collection
# end def mock_contigs_1

