'''
camprot.py - Tools for Computational Proteomics
===============================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Proteomics

To get help on a specific tool, type:

    camprot <tool> --help

To use a specific tool, type::

    camprot <tool> [tool options] [tool arguments]
'''

import os
import sys
import imp


def main(argv=None):

    argv = sys.argv

    path = os.path.abspath(os.path.dirname(__file__))

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print(globals()["__doc__"])

        return

    command = argv[1]

    (file, pathname, description) = imp.find_module(command, [path, ])
    module = imp.load_module(command, file, pathname, description)
    # remove 'umi-tools' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
