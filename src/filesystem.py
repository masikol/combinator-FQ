# -*- coding: utf-8 -*-

import re
import os

from src.platform import platf_depend_exit


def is_fasta(fpath: str) -> bool:
    # Returns True if path passed to it seems to point to a fasta file, else False.
    return not re.search(r'f(ast)?a(\.gz)?$', fpath) is None
# end def is_fasta


def make_outdir(outdpath: str) -> None:
    # Function creates output directory.
    if not os.path.exists(outdpath):
        try:
            os.makedirs(outdpath)
        except OSError as err:
            print('Error: cannot create output directory `{}`.'.format(outdpath))
            print(str(err))
            platf_depend_exit(1)
        # end try
    # end if
# end def make_outdir


def conf_prefix(infpath: str, outdpath: str) -> str:
    # Function returns output prefix for given input fasta file.

    # Make basic extention (without any numbers) by removing extention from file's name
    prefix: str = _rm_fasta_extention(infpath)

    # Check if output file corresponding to created prefix already exists.
    output_exists: bool = os.path.exists(
        os.path.join(outdpath, '{}_combinator_adjacent_contigs.tsv'.format(prefix))
    )

    # It it exists, we need to create new prefix by adding some number to basic one.
    if output_exists:
        # Pattern for "prefix with number"
        prefix_pattern: str = r'^{}\.([0-9]+)_combinator_adjacent_contigs\.tsv$'.format(prefix)
        # Count files with "prefix with number"
        num_prefix_files = len(tuple(
            filter(
                lambda f: not re.match(prefix_pattern, f) is None,
                os.listdir(outdpath)
            )
        ))

        # Guess that new prefix will have number `num_prefix_files + 1`
        curr_prefix_num: int = num_prefix_files + 1

        # Well, again, this file may exist
        output_exists = os.path.exists(
            os.path.join(outdpath, '{}_combinator_adjacent_contigs.tsv'\
                .format(_make_numbered_prefix(infpath, curr_prefix_num)))
        )

        # Increase "number" untill reaching non-extant output file
        while output_exists:
            curr_prefix_num += 1

            output_exists = os.path.exists(
                os.path.join(outdpath, '{}_combinator_adjacent_contigs.tsv'\
                    .format(_make_numbered_prefix(infpath, curr_prefix_num)))
            )

        # end while

        prefix = _make_numbered_prefix(infpath, curr_prefix_num)
    # end if

    return prefix
# end def make_prefix_if_exists


def _make_numbered_prefix(infpath: str, number: int) -> str:
    # Function makes "prefix with number" for given input file.
    return '{}.{}'.format(_rm_fasta_extention(infpath), number)
# end def _make_numbered_prefix


def _rm_fasta_extention(fpath: str) -> str:
    # Function removes fasta extention (with `.gz` one, if it it present)

    # Find the extention
    ext_match_obj: re.Match = re.search(
        r'.+(\.f(ast)?a(\.gz)?)$',
        os.path.basename(fpath)
    )

    # Remove it
    bname_no_ext: str
    if ext_match_obj is None:
        print('Error 12: please, contact the developer.')
        platf_depend_exit(12)
    else:
        bname_no_ext = os.path.basename(fpath).replace(ext_match_obj.group(1), '')
    # end if

    return bname_no_ext
# end def _rm_fasta_extention
