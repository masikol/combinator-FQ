# -*- coding: utf-8 -*-

import sys


def print_help(version: str, last_update_date: str) -> None:
    # Function prints help message.
    #
    # :param version: version of the program;
    # :param last_update_date: latest release last_update_date;

    print('    |=== {} ===|'.format(sys.argv[0]))
    print('Version {}. {} edition.\n'.format(version, last_update_date))
    print('Script identifies adjacent contigs in order to facilitate scaffolding.')
    print('Format of input: multi-fasta file containing contigs.')
    print('For assemblies made by SPAdes, combinator-FQ also considers coverage of each contig.')
    print("""\nScript generates 3 output files:
  1. Table containing table of adjecency.
Naming scheme: `<prefix>_combinator_adjacent_contigs.tsv`;
  2. File, in which all matches (not only adjacency-associated) are listed.
Naming scheme: `<prefix>_combinator_full_matching_log.txt`;
  3. Brief summary.
Naming scheme: `<prefix>_combinator_summary.txt`;""")
    print("""\n  Details can be found here:
https://github.com/masikol/combinator-FQ""")
    print('='*15 + '\n' + 'Options:\n')
    print('Length of an overlap is further referred to as `k`.\n')
    print('  -h (--help): print help message and exit.\n')
    print('  -v (--version): print version and exit.\n')
    print("""  -i (--mink): minimum k (in bp) to consider.
    Value: integer > 0; Default is 21 bp.\n""")
    print("""  -a (--maxk): maximum k (in bp) to consider.
    Value: integer > 0; Default is 127 bp.\n""")
    print("""  -k (--k-mer): exact k (in bp).
    If specified, `-i` and `-a` options are ignored.
    Value: integer > 0; Disabled by default.\n""")
    print("""  -o (--outdir): output directory.
    Default value: `combinator-result`.""")
    print('='*15 + '\n' + 'Examples:\n')
    print('  ./combinator-FQ.py contigs.fasta -k 127\n')
    print('  ./combinator-FQ.py another_contigs.fa -i 25 -a 300 -o my-outdir')
    print('\n'+'='*15)
    print("""If input file is omitted in the command, combinator-FQ will
  process all fasta files in the working directory.\n""")
# end def print_help
