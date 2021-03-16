# -*- coding: utf-8 -*-

import re


def is_fasta(fpath: str) -> bool:
    return not re.search(r'f(ast)?a(\.gz)?$', fpath) is None
# end def is_fasta
