#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import src.contigs
import src.overlaps
import src.statistics
import src.output

mink = 15
maxk = 25

contig_collection = src.contigs.get_contig_collection('test-contigs.fasta', maxk)
# contig_collection = src.contigs.get_contig_collection('/home/deynonih/Dropbox/tmp/contigs1000.fasta.gz', mink, maxk)
src.contigs.assign_multiplty(contig_collection)
overlap_collection = src.overlaps.detect_adjacent_contigs(contig_collection, mink, maxk)

# print(contig_collection)
print(overlap_collection[0])
print(overlap_collection[1])
print(overlap_collection[2])
print(overlap_collection[3])
print()

# print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 0))
# print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 1))
# print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 2))
# print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 3))
# print()

# print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 0))
# print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 1))
# print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 2))
# print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 3))
print()

print(src.statistics.calc_sum_contig_lengths(contig_collection))
cov_calc = src.statistics.CoverageCalculator(contig_collection)
print(cov_calc.get_min_coverage())
print(cov_calc.get_max_coverage())
print(cov_calc.calc_mean_coverage())
print(src.statistics.calc_lq_coef(contig_collection, overlap_collection))


src.output.write_summary(contig_collection, overlap_collection, 'infile', os.getcwd(), '')
src.output.write_adjacency_table(contig_collection, overlap_collection, os.getcwd(), '')
