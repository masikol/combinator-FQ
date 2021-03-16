#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import src.contigs
import src.overlaps
import src.statistics

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

print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 0))
print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 1))
print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 2))
print(src.overlaps.get_overlaps_str_for_table(overlap_collection, contig_collection, 3))
print()

print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 0))
print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 1))
print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 2))
print(src.overlaps.get_overlaps_str_for_log(overlap_collection, contig_collection, 3))
print()

print(src.statistics.calc_sum_contig_lengths(contig_collection))
print(src.statistics.get_min_coverage(contig_collection))
print(src.statistics.get_max_coverage(contig_collection))
print(src.statistics.calc_mean_coverage(contig_collection))
print(src.statistics.calc_lq_coef(contig_collection, overlap_collection))

