# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import getopt
from typing import List, Sequence, Dict, Mapping, Any, Tuple

import src.filesystem
from src.print_help import print_help
from src.platform import platf_depend_exit


def parse_args(version: str, last_update_date: str) -> Tuple[Sequence[str], Mapping[str, Any]]:

    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print_help(version, last_update_date)
        platf_depend_exit()
    # end if

    if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
        print(version)
        platf_depend_exit()
    # end if

    opts: List[List[str]]
    args: List[str]
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvk:i:a:o:',
            ['help', 'version', 'k-mer=', 'mink=', 'maxk=', 'outdir='])
    except getopt.GetoptError as err:
        print(str(err))
        platf_depend_exit(2)
    # end try

    contigs_fpaths: Sequence[str] = _get_input_fpaths(args)
    params: Dict[str, Any] = _parse_options(opts)

    # Verify mink and maxk:
    if params['i'] > params['a']:
        if not '-i' in sys.argv[1:] or not '--mink' in sys.argv[1:]:
            params['i'] = params['a']
        elif not '-a' in sys.argv[1:] or not '--maxk' in sys.argv[1:]:
            params['a'] = params['i']
        else:
            print('Error: minimum length of a k-mer is greater than maximum length of a k-mer.')
            print('Values specified by you:')
            print('Minimum length of a k-mer: {}.'.format(params['i']))
            print('Maximum length of a k-mer: {}.'.format(params['a']))
            platf_depend_exit(1)
        # end if
    # end if

    return contigs_fpaths, params
# end def parse_args


def _get_input_fpaths(args: Sequence[str]) -> Sequence[str]:

    contigs_fpaths: Sequence[str]

    # Determine fasta file to process:
    if len(args) == 0:
        # If no input file is specified,
        #   search for input files in the working directory:

        contigs_fpaths = tuple(
            filter(
                src.filesystem.is_fasta,
                glob.iglob(os.path.join(os.getcwd(), '*'))
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

    # Change paths to absolute paths
    contigs_fpaths = tuple(
        map(
            os.path.anspath,
            contigs_fpaths
        )
    )

    return contigs_fpaths
# end def _get_input_fpaths


def _parse_options(opts: Sequence[Sequence[str  ]]) -> Mapping[str, Any]:

    # Set default values for parameters
    params: Dict[str, Any] = {
        'o': os.path.join(os.getcwd(), 'combinator-result'), # outdir
        'i': 21,  # mink
        'a': 127, # maxk
    }

    # Parse command-line options

    for opt, arg in opts:

        if opt in ('-k', '--k-mer-len'):

            try:
                arg = int(arg)
                if arg <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: length of a k-mer must be positive integer number.')
                print('Your value: `{}`'.format(arg))
                platf_depend_exit(1)
            else:
                params['i'], params['a'] = arg, arg
            # end try

        elif opt in ('-o', '--outdir'):
            params['o'] = arg

        elif opt in ('-i', '--mink'):

            try:
                params['i'] = int(arg)
                if params['i'] <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Minimum length of a k-mer must be positive integer number.')
                print('Your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try

        elif opt in ('-a', '--maxk'):

            try:
                params['a'] = int(arg)
                if params['a'] <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Maximum length of k-mer must be positive integer number.')
                print('Your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try
        # end if
    # end for

    return params
# end def _parse_options
