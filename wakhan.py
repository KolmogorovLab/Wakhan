#!/usr/bin/env python3

#(c) 2025 by Authors
#This file is a part of Wakhan program.
#Released under the MIT license (see LICENSE file)

"""
This script sets up environment paths
and invokes Wakhan without installation.
"""

import os
import sys

def main():
    #Setting executable paths
    wakhan_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, wakhan_root)
    #bin_absolute = os.path.join(severus_root, BIN_DIR)
    #os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]

    #Severus entry point
    from src.main import main
    sys.exit(main())


if __name__ == "__main__":
    main()