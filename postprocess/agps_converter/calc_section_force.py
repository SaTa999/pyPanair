#!/usr/bin/env python
__author__ = "stakanashi"
import numpy as np
import pandas as pd
from .agps_reader import read_agps


def rot(v, theta, degrees=True):
    if degrees:
        theta = np.radians(theta)
    r = np.array([[np.cos(-theta), -np.sin(-theta)],
                  [np.sin(-theta), np.cos(-theta)]])
    return np.dot(r, v)


def calc_section_force(aoa, mac, centerofgravity, casenum, networknum):
    """calculate local aerodynamic coefficients from the pressure distribution in agps files
    aoa: AoA of the case (e.g. aoa = 2.)
    mac: mean aerodynamic chord of the wing

    reference: https://cfdws.chofu.jaxa.jp/apc/"""
    net_num = 1
    coordinates = ["x", "z"]
    span = "y"
    coor2id = {"x":1, "y":2, "z":3}
    # 迎角、翼弦長等を読み込む

    rot_center = np.asarray(centerofgravity)

    # read the pressure distribution from agps
    dat = read_agps()
    network = dat[net_num-1]

    result = []

    for data in network:
        data = data[:-1]
        # 断面圧力分布を読み込み、局所空力係数を返す
        y = float(data[0,coor2id[span]])
        data_clip = data[:, [coor2id[c] for c in coordinates]]
        diff_coord = (data_clip - np.roll(data_clip, 1, axis=0))
        l = np.linalg.norm(diff_coord, axis=1) # 格子点間の距離
        norm = np.fliplr(diff_coord / l[:,np.newaxis]) # 法線ベクトル
        min_id = np.argmin(data_clip, axis=0) # 座標が最小となる行
        flip = np.array([-1 if col[id] < 0 else 1 for (col, id) in zip(norm.T, min_id)])
        norm = norm * flip # 法線ベクトルを内側に向ける
        chord = np.max(data_clip[:,0]) - np.min(data_clip[:,0])

        coeff_all = [y]
        for (cp, aoa) in zip(data.T[4:], aoas):
            cp2 = cp
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
    columns = ["pos"]
    columns += "cd cl cm".split()
    result = pd.DataFrame(result, columns=columns)
    result.to_csv("section_force.csv", index=False)