#!/usr/bin/env python
__author__ = "stakanashi"
import numpy as np
import pandas as pd
from pyPanair.postprocess import read_agps


def rot(vec, angle, degrees=True):
    """rotate the vector for "angle" degrees (or radians)"""
    if degrees:
        angle = np.radians(angle)
    r = np.array([[np.cos(-angle), -np.sin(-angle)],
                  [np.sin(-angle), np.cos(-angle)]])
    return r @ vec


def calc_section_force(aoa, mac, rot_center, casenum=1, networknum=1):
    """calculate local aerodynamic coefficients from the pressure distribution in agps files
    :param aoa: AoA of the case in degrees (e.g. 2.94)
    :param mac: mean aerodynamic chord of the wing (e.g. 275.8 [in])
    :param rot_center: center of rotation of the aircraft (e.g. (4.80, 0.00, 0.65) [in])
    :param casenum: case number to calculate the local aerodynamic coefficients (e.g. 1)
    :param networknum: the network number of the wing (e.g. 1)
    reference: JAXA, "APC2 Submission Format," https://cfdws.chofu.jaxa.jp/apc/apc2/submit.html"""

    chord = "x"
    span = "y"
    height = "z"
    coor2id = {"x":1, "y":2, "z":3}
    chord_id = coor2id[chord]
    span_id = coor2id[span]
    height_id = coor2id[height]

    rot_center = np.asarray(rot_center)[[chord_id-1, height_id-1]]

    # read the pressure distribution from agps
    data = read_agps()
    network = data[networknum-1]
    result = []

    # calculate the local aerodynamic coefficients for each line in the network
    for line in network:
        line = line[:-1] # omit the last vertex in the line (because it has the same coordinate as the first vertex)
        y = float(line[0,span_id]) # coordinate of spanwise position
        line_clip = line[:, [chord_id, height_id]] # coordinate of cross section
        diff_coord = line_clip - np.roll(line_clip, 1, axis=0)
        length = np.linalg.norm(diff_coord, axis=1) # length between each vertex
        norm = np.fliplr(diff_coord / length[:,np.newaxis]) # norm vector between each vertex
        min_id = np.argmin(line_clip, axis=0) # vertex_id with the minimum coordinate
        flip = np.array([-1 if col[id] < 0 else 1 for (col, id) in zip(norm.T, min_id)])
        norm *= flip # flip the norm vectors to make them point inward
        chord = np.linalg.norm(line_clip[0] - line_clip[line_clip.shape[0]//2])

        coeff_all = [y]
        # the definition of each variable is explained in the reference
        cp2 = line.T[3+casenum]
        cp1 = np.roll(cp2, 1, axis=0)
        cp3 = np.roll(cp2, -1, axis=0)
        n21 = norm
        n23 = np.roll(n21, -1, axis=0)
        l21 = length * 0.5
        l23 = np.roll(l21, -1, axis=0)
        cp21 = (0.75 * cp2 + 0.25 * cp1)[:,np.newaxis] * n21
        cp23 = (0.75 * cp2 + 0.25 * cp3)[:,np.newaxis] * n23
        cpdl = cp21 * l21[:,np.newaxis] + cp23 * l23[:,np.newaxis]
        cp_moment = cpdl * np.fliplr(line_clip - rot_center)
        cp_moment = cp_moment[:,0] - cp_moment[:,1]
        liftdragcoef = np.sum(cpdl, axis=0) / mac
        liftdragcoef = rot(liftdragcoef, aoa)
        cm = np.sum(cp_moment, axis=0) / (mac * chord)
        coeff = liftdragcoef.tolist()
        coeff.append(cm)
        coeff_all += coeff
        result.append(coeff_all)
    result = np.array(result)
    columns = ["pos", "cd", "cl", "cm"]
    result = pd.DataFrame(result, columns=columns)
    result.to_csv("section_force.csv", index=False)