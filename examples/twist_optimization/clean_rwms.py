#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
from os import remove
from time import sleep


def main(target_dir):
    # rwmsを削除
    rwms = glob("{}/rwms*".format(target_dir))
    while rwms:
        for f in rwms:
            try:
                remove(f)
            except WindowsError:
                print("failed to delete {}".format(f))
        rwms = glob("{}/rwms*".format(target_dir))
        if rwms:
            sleep(2)


