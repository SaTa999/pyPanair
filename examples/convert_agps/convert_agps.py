from pyPanair.postprocess import agps_converter


def inputloop(message, outtype):
    while True:
        try:
            print(message, end="")
            inp = outtype(input())
            return inp
        except ValueError:
            print("invalid input, try again")


if __name__ == '__main__':
    # test conversion to vtk, vtm, and dat
    print("select the format of the output file")
    print("1: vtk (legacy paraview format)")
    print("2: vtm (multiblock paraview format)")
    print("3: dat (mutliblock tecplot format)")
    while True:
        out_format = inputloop("enter 1, 2, or 3: ", int)
        if 1 <= out_format <= 3:
            break
        else:
            print("invalid input, try again")
    print("==============================================")
    print("specify the number of wakes")
    print("the last \"n_wakes\" networks in the agps file  will not be included in the output file")
    while True:
        n_wakes = inputloop("enter a number: ", int)
        if n_wakes >= 0:
            break
        else:
            print("invalid input, try again")
    print("==============================================")
    switch_dict = {1: agps_converter.write_vtk, 2: agps_converter.write_vtm, 3: agps_converter.write_tec}
    switch_dict[out_format](n_wakes)
