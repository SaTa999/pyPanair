#!/usr/bin/env python
import os
import numpy as np
from scipy.interpolate import interp1d
from pyPanair.preprocess import wgs_creator


def main(twist_func, aoas=(7.42), target_dir=""):
    """ create a LaWGS file for twisted rectangular wing
    reference case 3"""
    wgs = wgs_creator.LaWGS("ADODG_case3")
    n_wing_x = 40
    n_wing_y = 40
    n_wing_z = 12
    chord_wing = 100.
    halfspan_wing = chord_wing * 3

    # define twist distribution
    twist_xy = twist_func(np.linspace(0, 1, 150))
    twist_dist = interp1d(twist_xy[:,0] * halfspan_wing, twist_xy[:,1])

    # create wing
    base_airfoil = wgs_creator.naca4digit("0010", num=n_wing_x, chord=chord_wing)
    wing = list()
    span_pos = np.linspace(0, halfspan_wing, n_wing_y)
    for y in span_pos:
        rot_foil = base_airfoil.roty(base_airfoil[0], twist_dist(y))
        rot_foil.shift((0, y, 0), inplace=True)
        wing.append(rot_foil)
    wing = wgs_creator.Network(wing)
    wgs.append_network("wing", wing, 1)

    # create wingtip
    degs = np.linspace(0, -180, n_wing_z)
    wingtipu, _ = base_airfoil.shift((0, halfspan_wing, 0)).split_half()
    wingtip = [wingtipu.rotx(wingtipu[0], d) for d in degs]
    wingtip = wgs_creator.Network(wingtip)
    wingtip = wingtip.roty(wingtipu[0], twist_dist(halfspan_wing))
    wgs.append_network("wingtip", wingtip, 1)

    # add wake
    wake_length = chord_wing * 50.
    wingwake = wing.make_wake(3, wake_length)
    wgs.append_network("wingwake", wingwake, 18)

    wgs.create_wgs(os.path.join(target_dir, "ADODG_case3.wgs"))

    span = halfspan_wing * 2
    sref = chord_wing * span / 100 # multiply aerodynamic coefficients by 100
    xref = chord_wing * 0.25
    aux_name = os.path.join(target_dir, "ADODG_case3.aux")
    wgs.create_aux(filename=aux_name, alpha=aoas, mach=0.5, cbar=chord_wing, span=span, sref=sref,
                   xref=xref, zref=0.)

    # wgs.create_stl("ADODG3.stl")


if __name__ == '__main__':
    main()
