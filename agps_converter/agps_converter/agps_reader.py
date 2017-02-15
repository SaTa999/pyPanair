#!/usr/bin/env python
__author__ = "stakanashi"
import numpy as np

def read_column(file, firstline):
    """read a column from first line (e.g. n01c001) to *eof"""
    column = list()
    line = firstline
    # register methods for faster evaluation
    f_readline = file.readline
    column_append = column.append
    # read each line until *eof
    while line:
        line = f_readline().split()
        if line[0] == "*eof":
            break
        column_append(line)
    return column


def read_network(file, firstheader):
    """read a network"""
    network_n = int(firstheader[0][1:3]) # get network number from the first header (e.g. 01 from n01c001)
    print("loading network no.", network_n)
    network = list()
    line = firstheader
    # register methods for faster evaluation
    network_append = network.append
    # read each line until next header
    while line:
        col = read_column(file, line)
        network.append(col)
        line = file.readline().split()
        # break at the end of agps file
        if not line:
            break
        # break when reaching the header for the next network (e.g. n02c001)
        if not int(line[0][1:3]) == network_n:
            break
    network = np.array(network, dtype=float)
    return network, line


def read_agps():
    # read the agps file and return a list of arrays containing data for each network
    with open("agps", "r") as f:
        # skip the header of the agps file
        for _ in range(6):
            line = f.readline()
        line = f.readline().split()
        f.readline() # skip the header of first network ("icol, x, y, z, cp1, cp2, cp3, cp4")
        dat = []
        while line:
            net, line = read_network(f, line)
            dat.append(net)
    return dat
