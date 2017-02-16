#!/usr/bin/env python
from .ffmf_reader import read_ffmf

# def natsort(l):
#   convert = lambda text: int(text) if text.isdigit() else text
#   alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
#   l.sort( key=alphanum_key )
  
def write_ffmf(outfilepath="ffmf.csv"):
    ffmf = read_ffmf()
    ffmf.to_csv(outfilepath, index=None)
