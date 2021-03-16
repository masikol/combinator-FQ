#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

__version__ = '1.4.a'
# Year, month, day
__last_update_date__ = '2021-03-XX'
__min_python_version__ = 3.6

# Check python interpreter version
if sys.version_info.major + sys.version_info.minor*0.1 < __min_python_version__:
    print( 'Your python interpreter version is ' + '%d.%d' % (sys.version_info.major,
        sys.version_info.minor) )
    print('  Please, use Python %.1f+.\a' % __min_python_version__)
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        if sys.version_info.major == 2:
            raw_input('Press ENTER to exit:')
        else:
            input('Press ENTER to exit:')
        # end if
    # end if
    sys.exit(1)
# end if

from src.main import main

if __name__ == '__main__':
    main()
# end if
