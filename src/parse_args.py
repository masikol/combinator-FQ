# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import getopt
from typing import List, Collection

import src.filesystem
from src.print_help import print_help
from src.platform import platf_depend_exit


def parse_args():

    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print_help()
        platf_depend_exit()
    # end if

    if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
        print(__version__)
        platf_depend_exit()
    # end if

    opts: List[List[str]]
    args: List[str]
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvp:k:i:a:o:',
            ['help', 'version', 'prefix=', 'k-mer-len=', 'mink=', 'maxk=', 'outdir='])
    except getopt.GetoptError as err:
        print( str(err) )
        platf_depend_exit(2)
    # end try

    contigs_fpaths = _get_input_fpaths(args)

    print(contigs_fpaths)

    # prefix = re.search(r"(.+)\.(m)?f(ast)?a(\.gz)?", os.path.basename(contigs_fpath)).group(1)
    # outdpath = "combinator-result"
    # mink = 21
    # maxk = 127

    # # Parse command-line options

    # for opt, arg in opts:

    #     if opt in ("-k", "--k-mer-len"):

    #         try:
    #             arg = int(arg)
    #             if arg <= 0:
    #                 raise ValueError
    #             # end if
    #         except ValueError:
    #             print("\nLength of k-mer must be positive integer number.")
    #             print("Your value: '{}'".format(arg))
    #             platf_depend_exit(1)
    #         else:
    #             mink, maxk = arg, arg # mink = k and maxk = k
    #         # end try

    #     elif opt in ("-o", "--outdir"):
    #         outdpath = arg

    #     elif opt in ("-i", "--mink"):

    #         try:
    #             mink = int(arg)
    #             if mink <= 0:
    #                 raise ValueError
    #             # end if
    #         except ValueError:
    #             print("\nMinimum length of k-mer must be positive integer number.")
    #             print("Your value: '{}'".format(arg))
    #             platf_depend_exit(1)
    #         # end try

    #     elif opt in ("-a", "--maxk"):

    #         try:
    #             maxk = int(arg)
    #             if maxk <= 0:
    #                 raise ValueError
    #             # end if
    #         except ValueError:
    #             print("\nMaximum length of k-mer must be positive integer number.")
    #             print("Your value: '{}'".format(arg))
    #             platf_depend_exit(1)
    #         # end try

    #     elif opt in ("-p", "--prefix"):
    #         prefix = arg
    #     # end if
    # # end for

    # # Verify mink and maxk:
    # if mink > maxk:
    #     if not "-i" in sys.argv[1:] or not "--mink" in sys.argv[1:]:
    #         mink = maxk
    #     elif not "-a" in sys.argv[1:] or not "--maxk" in sys.argv[1:]:
    #         maxk = mink
    #     else:
    #         print(err_fmt("Minimum length of k-mer is greater than maximum length of k-mer."))
    #         print("Values specified by you:")
    #         print("Minimum length of k-mer: {}.".format(mink))
    #         print("Maximum length of k-mer: {}.".format(maxk))
    #         platf_depend_exit(1)
    #     # end if
    # # end if
# end def parse_args


def _get_input_fpaths(args: Collection[str]) -> Collection[str]:

    contigs_fpaths: Collection[str]

    # Determine fasta file to process:
    if len(args) == 0:
        # If no input file is specified,
        #   search for input files in the working directory:

        contigs_fpaths = tuple(
            filter(
                src.filesystem.is_fasta,
                glob.iglob(os.path.join(os.getcwd(), '*contigs.f*'))
            )
        )

        if len(contigs_fpaths) != 0:
            pass
        else:
            # If there are no fasta files in working directory, just print help
            print('Nothing to process.')
            print('Please, type `{} -h` for help.'.format(sys.argv[0]))
            platf_depend_exit(1)
        # end if

    else:
        # Check existance of input files
        fpath: str
        for arg in args:
            if not os.path.exists(arg):
                print('Error: file `{}` does not exist.'.format(arg))
                platf_depend_exit(1)
            # end if
            if not src.filesystem.is_fasta(arg):
                print('Error: file `{}` is not a fasta file, considering it\'s extention.'.format(arg))
                platf_depend_exit(1)
            # end if
        # end for
        contigs_fpaths = tuple(args)
    # end if

    return contigs_fpaths
# end def _get_input_fpaths
