import agps_converter

if __name__ == '__main__':
    # test conversion to vtk, vtm, and dat
    agps_converter.write_vtk(0)
    agps_converter.write_vtm(0)
    agps_converter.write_tec(0)
