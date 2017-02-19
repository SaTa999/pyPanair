from pyPanair.postprocess import calc_section_force


def inputloop(message, outtype):
    while True:
        try:
            print(message, end="")
            inp = outtype(input())
            return inp
        except ValueError:
            print("invalid input, try again")


if __name__ == '__main__':
    # test calculation of section forces
    aoa = inputloop("enter the AoA of the case: ", float)
    mac = inputloop("enter the mac of the wing: ", float)
    print("enter the center of rotation")
    rot_center = ()
    for cor in ("x", "y", "z"):
        rot_center += (inputloop("{}: ".format(cor), float), )
    casenum = inputloop("enter the case number: ", int)
    networknum = inputloop("enter the network number of the wing: ", int)

    calc_section_force(aoa, mac, rot_center, casenum, networknum)
