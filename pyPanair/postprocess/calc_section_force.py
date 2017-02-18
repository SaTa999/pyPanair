#!/usr/bin/env python
__author__ = "stakanashi"
import numpy as np
import pandas as pd
from .agps_converter import read_agps


def rot(v, theta, degrees=True):
    """rotate the vector for theta degrees (or radians)"""
    if degrees:
        theta = np.radians(theta)
    r = np.array([[np.cos(-theta), -np.sin(-theta)],
                  [np.sin(-theta), np.cos(-theta)]])
    return r @ v


def calc_section_force(aoa, mac, rot_center, casenum, networknum):
    """calculate local aerodynamic coefficients from the pressure distribution in agps files
    :param aoa: AoA of the case in degrees (e.g. 2.94)
    :param mac: mean aerodynamic chord of the wing (e.g. 275.8 [in])
    :param rot_center: center of rotation of the aircraft (e.g. (4.80, 0.00, 0.65) [in])
    :param casenum: case number to calculate the local aerodynamic coefficients (e.g. 1)
    :param networknum: the network number of the wing (e.g. 1)
    reference: JAXA, "APC2 Submission Format," https://cfdws.chofu.jaxa.jp/apc/apc2/submit.html"""
    rot_center = np.asarray(rot_center)

    coordinates = ["x", "z"]
    span = "y"
    coor2id = {"x":1, "y":2, "z":3}

    # read the pressure distribution from agps
    dat = read_agps()
    network = dat[networknum-1]
    result = []

    # calculate the local aerodynamic coefficients for each line in the network
    for line in network:
        line = line[:-1]
        y = float(line[0,coor2id[span]])
        line_clip = line[:, [coor2id[c] for c in coordinates]]
        diff_coord = (line_clip - np.roll(line_clip, 1, axis=0))
        l = np.linalg.norm(diff_coord, axis=1) # 格子点間の距離
        norm = np.fliplr(diff_coord / l[:,np.newaxis]) # 法線ベクトル
        min_id = np.argmin(line_clip, axis=0) # 座標が最小となる行
        flip = np.array([-1 if col[id] < 0 else 1 for (col, id) in zip(norm.T, min_id)])
        norm = norm * flip # 法線ベクトルを内側に向ける
        chord = np.max(line_clip[:,0]) - np.min(line_clip[:,0])

        coeff_all = [y]
        cp2 = data.T[3+casenum:]
        cp1 = np.roll(cp2, 1, axis=0)
        cp3 = np.roll(cp2, -1, axis=0)
        n21 = norm
        n23 = np.roll(n21, -1, axis=0)
        l21 = l * 0.5
        l23 = np.roll(l, -1, axis=0) * 0.5
        cp21 = (0.75 * cp2 + 0.25 * cp1)[:,np.newaxis] * n21
        cp23 = (0.75 * cp2 + 0.25 * cp3)[:,np.newaxis] * n23
        cpdl = cp21 * l21[:,np.newaxis] + cp23 * l23[:,np.newaxis]
        cp_moment = cpdl * np.fliplr(data_clip - rot_center)
        cp_moment = cp_moment[:,0] - cp_moment[:,1]
        coeff = np.sum(cpdl, axis=0) / mac
        coeff = rot(coeff, aoa)
        cm = np.sum(cp_moment, axis=0) / mac / chord
        coeff = coeff.tolist()
        coeff.append(cm)
        coeff_all += coeff
        result.append(coeff_all)
    result = np.array(result)
    columns = ["pos", "cd", "cl", "cm"]
    result = pd.DataFrame(result, columns=columns)
    result.to_csv("section_force.csv", index=False)