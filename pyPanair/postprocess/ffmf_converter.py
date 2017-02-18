#!/usr/bin/env python
__author__ = "stakanashi"
import pandas as pd


def read_ffmf(inpfilepath="ffmf"):
    ffmf = []
    columns = ["sol-no", "AoA", "beta", "cl", "cdi", "cy", "fx", "fy", "fz", "mx", "my", "mz", "area"]
    sol_num = 1
    with open(inpfilepath, "r") as f:
        for i in range(30): # force / torque coefficients are listed in the first 30 rows
            line = f.readline().split()
            if line == []:
                pass
            # read force / torque coefficients
            elif line[0] == str(sol_num):
                sol_num += 1
                line2 = f.readline().split()
                ffmf.append(line + line2)
            else:
                pass
    ffmf = pd.DataFrame(ffmf, columns=columns)
    return ffmf

  
def write_ffmf(outfilepath="ffmf.csv"):
    ffmf = read_ffmf()
    ffmf.to_csv(outfilepath, index=None)
