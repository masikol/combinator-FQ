# -*- coding: utf-8 -*-

import os
import sys
import glob
import getopt
from typing import List, Sequence, Dict, Mapping, Any, Tuple

import src.filesystem
from src.print_help import print_help
from src.platform import platf_depend_exit


def parse_args(version: str, last_update_date: str) -> Tuple[Sequence[str], Mapping[str, Any]]:
    # Function parses command line arguments.
    # Returns two values:
    #  1. Collection of paths to input files.
    #  2. Dictionary of parameters (see function _parse_options).

    # Print help message and exit if required
    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print_help(version, last_update_date)
        platf_depend_exit()
    # end if

    # Print version and exit if required
    if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
        print(version)
        platf_depend_exit()
    # end if

    # Parse arguments woth getopt
    opts: List[List[str]]
    args: List[str]
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvk:i:a:o:',
            ['help', 'version', 'k-mer=', 'mink=', 'maxk=', 'outdir='])
    except getopt.GetoptError as err:
        print(str(err))
        platf_depend_exit(2)
    # end try

    # Extract paths to input files from parsed arguments
    contigs_fpaths: Sequence[str] = _get_input_fpaths(args)
    # Extract optional parameters from parsed arguments
    params: Dict[str, Any] = _parse_options(opts)

    # Verify mink and maxk:
    if params['i'] > params['a']:
        if '-i' not in sys.argv[1:] or '--mink' not in sys.argv[1:]:
            params['i'] = params['a']
        elif '-a' not in sys.argv[1:] or '--maxk' not in sys.argv[1:]:
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
    # Function extracts paths to input files from `args` colection returned by `getopt.gnu_getopt`.
    # Returns collection of paths to input fasta files.

    contigs_fpaths: Sequence[str]

    if len(args) == 0:
        # If no input files are specified,
        #   search for input files in the working directory:

        # Get all fasta files in the workgin directory
        contigs_fpaths = tuple(
            filter(
                src.filesystem.is_fasta,
                glob.iglob(os.path.join(os.getcwd(), '*'))
            )
        )

        if len(contigs_fpaths) != 0:
            # Files are found -- ask for permission and proceed
            print('Following fasta files are found and will be processed:')
            i: int
            path: str
            for i, path in enumerate(contigs_fpaths):
                print(' {}. `{}`'.format(i+1, path))
            # end for
            error: bool = True
            while error:
                reply: str = input('Is it ok? Proceed? [y, n]:')
                if reply.lower() == 'y':
                    error = False
                elif reply.lower() == 'n':
                    print('Ok. Type `{} -h` for help.'.format(sys.argv[0]))
                    sys.exit(0)
                else:
                    print('Invalid reply: `{}`'.format(reply))
                # end if
            # end while
        else:
            # If there are no fasta files in the working directory -- exit
            print('Nothing to process.')
            print('Please, type `{} -h` for help.'.format(sys.argv[0]))
            platf_depend_exit(1)
        # end if
    else:
        # Check existance of input files
        arg: str
        for arg in args:
            if not os.path.exists(arg):
                print('Error: file `{}` does not exist.'.format(arg))
                platf_depend_exit(1)
            # end if
            if not src.filesystem.is_fasta(arg):
                print('Error: file `{}` is not a fasta file, considering it\'s extention.'\
                    .format(arg))
                platf_depend_exit(1)
            # end if
        # end for
        contigs_fpaths = tuple(args)
    # end if

    # make all paths absolute
    contigs_fpaths = tuple(
        map(
            os.path.abspath,
            contigs_fpaths
        )
    )

    return contigs_fpaths
# end def _get_input_fpaths


def _parse_options(opts: Sequence[Sequence[str]]) -> Mapping[str, Any]:
    # Funtion extracts options from `args` colection returned by `getopt.gnu_getopt`.
    # Returns following dictionary:
    #    {
    #       'o': <outdir_path>,
    #       'i': <mink>,
    #       'a': <maxk>,
    #    }

    # Set default values for parameters
    params: Dict[str, Any] = {
        'o': os.path.join(os.getcwd(), 'combinator-result'), # outdir
        'i': 21,                                             # mink
        'a': 127,                                            # maxk
    }

    # Parse command line options
    opt: str
    arg: str
    for opt, arg in opts:

        # k-mer size
        if opt in ('-k', '--k-mer-len'):
            try:
                k: int = int(arg)
                if k <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: length of a k-mer must be positive integer number.')
                print('Your value: `{}`'.format(arg))
                platf_depend_exit(1)
            else:
                params['i'], params['a'] = k, k
            # end try

        # Outdir
        elif opt in ('-o', '--outdir'):
            params['o'] = arg

        # Minimim k-mer size
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

        # Maximim k-mer size
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
