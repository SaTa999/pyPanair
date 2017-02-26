#!/usr/bin/env python
import glob
import os
import re
from pyPanair.postprocess import write_vtk


def convert(text):
    return int(text) if text.isdigit() else text


def alphanum_key(key):
    return [convert(c) for c in re.split('([0-9]+)', key)]


def natsort(l):
    l.sort(key=alphanum_key)


if __name__ == '__main__':
    results = glob.glob("results/case*/agps")
    natsort(results)
    if not os.path.exists("video"):
        os.mkdir("video")
    re_num = re.compile(r"case\d+")
    for r in results:
        case = re_num.search(r).group(0)
        write_vtk(n_wake=1, inputfile=r, outputname=os.path.join("video", case))
