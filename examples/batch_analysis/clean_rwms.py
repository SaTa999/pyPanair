#!/usr/bin/env python
from glob import glob
from os import remove
from warnings import warn


if __name__ == '__main__':
    rwms = glob("rwms*")
    while rwms:
        for f in rwms:
            try:
                remove(f)
            except OSError:
                warn("failed to delete {}".format(f))
        rwms = glob("rwms*")

