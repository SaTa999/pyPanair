#!/usr/bin/env python
from functools import partial
import numpy as np
from scipy.interpolate import interp1d
from pyPanair.preprocess import wgs_creator
from pyPanair.utilities import bspline


def main(target_dir=None, twist_func=None, aoas=(7.42)):
    """ create a LaWGS file for twisted rectangular wing
    reference case 3"""
    wgs = wgs_creator.LaWGS("model for ADODG case3")
    n_wing_x = 40
    n_wing_y = 40
    n_wing_z = 12
    chord_wing = 1.
    halfspan_wing = chord_wing * 3

    cv = np.array([[0., 0.],
                   [0., -4.999691458172163649],
                   [0.5000349654015733281, -2.231521636034730971],
                   [1., -4.732433281976027750]])

    # n_cpoints = 4
    # sep = 1 / (n_cpoints - 2)
    # cv = np.zeros((n_cpoints, 2))
    # for i in range(1, n_cpoints-1):
    #     cv[i,0] = np.random.random() * sep + sep * (i-1)
    # cv[-1,0] = 1
    # cv[1:,1] = np.random.random(cv.shape[0]-1) * 20 -10
    twist_func = bspline(cv, degree=3, periodic=False)

    twist_xy = twist_func(np.linspace(0, 1, 150))
    twist_dist = interp1d(twist_xy[:,0] * halfspan_wing, twist_xy[:,1])

    base_airfoil = wgs_creator.naca4digit("0010", num=n_wing_x, chord=chord_wing)
    wing = list()
    span_pos = np.linspace(0, halfspan_wing, n_wing_y)
    for y in span_pos:
        rot_foil = base_airfoil.roty(base_airfoil[0], twist_dist(y))
        rot_foil.shift((0, y, 0), inplace=True)
        wing.append(rot_foil)
    wingtipu, _ = wing[-1].split_half()
    wing = wgs_creator.Network(wing)
    wgs.append_network("wing", wing, 1)

    # wingtip = wingtipu.linspace(wingtipl.flip(), n_wing_z) # flat tip
    degs = np.linspace(0, -180, n_wing_z)
    wingtip = [wingtipu.rotx(wingtipu[0], d) for d in degs]
    wingtip = wgs_creator.Network(wingtip)
    wingtip = wingtip.roty(wingtipu[0], twist_dist(halfspan_wing))
    wgs.append_network("wingtip", wingtip, 1)

    # add wake
    wake_length = chord_wing * 50.
    wingwake = wing.make_wake(3, wake_length)
    wgs.append_network("wingwake", wingwake, 18)

    wgs.create_wgs("ADODG3.wgs".format(target_dir))

    span = halfspan_wing * 2
    sref = chord_wing * span
    xref = chord_wing * 0.25
    wgs.create_aux("ADODG3.aux".format(target_dir), alpha=aoas, mach=0.5, cbar=chord_wing, span=span, sref=sref,
                   xref=xref, zref=0.)

    wgs.create_stl("ADODG3.stl")


if __name__ == '__main__':
    main()
