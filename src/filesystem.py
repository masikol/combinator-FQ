# -*- coding: utf-8 -*-

import re
import os

from src.platform import platf_depend_exit


def is_fasta(fpath: str) -> bool:
    return not re.search(r'f(ast)?a(\.gz)?$', fpath) is None
# end def is_fasta


def make_outdir(outdpath: str) -> None:
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
