# -*- coding: utf-8 -*-

import sys
from typing import NoReturn

def platf_depend_exit(exit_code: int = 0) -> NoReturn:
    # Function asks to press ENTER press on Windows
    #     and exits after that.
    # :type exit_code: int;

    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit
