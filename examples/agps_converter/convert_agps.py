from pyPanair.postprocess import agps_converter

if __name__ == '__main__':
    # test conversion to vtk, vtm, and dat
    print("select the format of the output file")
    print("1: vtk (legacy paraview format)")
    print("2: vtm (multiblock paraview format)")
    print("3: dat (mutliblock tecplot format)")
    print("enter 1, 2, or 3: ", end="")
    while True:
        try:
            out_format = int(input())
            if 1 <= out_format <= 3:
                break
            else:
                print("enter 1, 2, or 3: ", end="")
        except ValueError:
            print("enter 1, 2, or 3: ", end="")

    print()
    print("specify the number of wakes")
    print("the last \"n_wakes\" networks in the agps file  will not be included in the output file")
    print("enter a number: ", end="")
    while True:
        try:
            n_wakes = int(input())
            break
        except ValueError:
            print("enter a number: ", end="")

    print("==============================================")
    switch_dict = {1: agps_converter.write_vtk, 2: agps_converter.write_vtm, 3: agps_converter.write_tec}
    switch_dict[out_format](n_wakes)
