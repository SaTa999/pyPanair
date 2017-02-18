#!/usr/bin/env python
__author__ = "stakanashi"
import numpy as np
import os


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


def read_agps(inputfile="agps"):
    # read the agps file and return a list of arrays containing data for each network
    with open(inputfile, "r") as f:
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


def write_vtk(n_wake=0, outputname="agps", inputfile="agps"):
    """Write agps in the legacy paraview format (vtk)
    All networks will be merged into one block
    Therefore, user are advised to omit 'wakes' by specifying the 'n_wakes'"""
    data = read_agps(inputfile) # read agps file & specify the number of networks to omit
    print("n_wake = ", n_wake)
    # write the header of the vtk file
    vtk = "# vtk DataFile Version 2.0\n"
    vtk += "scalar\n"
    vtk += "ASCII\n"
    vtk += "DATASET UNSTRUCTURED_GRID\n"
    n_points = 0 # number of points in vtk file
    n_cells = 0 # number of quadrilateral cells formed by the points
    n_cp = data[0].shape[2] - 4
    points = str() # coordinate of each point (x, y, z)
    point_data = [str()] * n_cp # cp at each point (cp1, cp2, cp3, cp4)
    cells = str() # ids of each quadrilateral cell (e.g. (0, n_col, n_col + 1, 1) for first cell)
    for i in range(len(data)-n_wake):
        net = data[i]
        n_row = int(net.shape[0])
        n_col = int(net.shape[1])
        print("network {} shape: ".format(i+1), net.shape)
        base_square = np.array((0, n_col, n_col+1, 1))
        for j in range(n_row):
            for k in range(n_col):
                point = net[j,k]
                # add coordinate of a point
                points += "{0} {1} {2}\n".format(point[1], point[2], point[3])
                # add cp data of a point
                for l in range(n_cp):
                    point_data[l] += "{}\n".format(point[4+l])
                # add ids of a cell
                if not j == n_row-1 and not k == n_col-1:
                    square =  base_square + (j * n_col + k) + n_points
                    square = (str(p) for p in square)
                    cells += "4 " + " ".join(square) + "\n"
        # add the number of points / cells
        n_points += n_row * n_col
        n_cells += (n_row - 1) * (n_col - 1)
    # write the header of each block (POINTS, CELLS, CELLTYPES, POINT_DATA)
    points = "POINTS {} float\n".format(n_points) + points
    cells = "CELLS {0} {1}\n".format(n_cells, n_cells*5) + cells
    cell_types = "CELL_TYPES {}\n".format(n_cells) + "9\n"*n_cells
    vtk += points + cells + cell_types + "POINT_DATA {}\n".format(n_points)
    for l in range(n_cp):
        vtk += "SCALARS cp{} float\nLOOKUP_TABLE default\n".format(l+1) + point_data[l]
    with open("{}.vtk".format(outputname), "w") as f:
        f.write(vtk)


def write_vtm(n_wake=0, outputname="agps", inputfile="agps"):
    """convert agps networks to paraview unstructured grid
    each network will become a different vtu file
    to open all vtu files at the same time, open the vtm file with paraview"""
    data = read_agps(inputfile) # read agps file & specify the number of networks to omit
    print("n_wake = ", n_wake)
    # write header of vtm file
    vtm = "<?xml version=\"1.0\"?>\n"
    vtm += "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
    vtm += "  <vtkMultiBlockDataSet>\n"

    for i in range(len(data)-n_wake):
        # add dataset to vtm file
        vtu_dir = "{}_vtu".format(outputname)
        try:
            os.mkdir(vtu_dir)
        except OSError:
            if not os.path.exists(vtu_dir):
                raise
        vtu_path = "{0}/{1}{2}.vtu".format(vtu_dir, outputname, i+1)
        vtm += "    <DataSet index=\"network{0}\" file=\"{1}\"/>\n".format(i+1, vtu_path)
        # write header of vtu file
        vtu = "<?xml version=\"1.0\"?>\n"
        vtu += "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        vtu += "  <UnstructuredGrid>\n"
        # write the header of the piece
        net = data[i]
        n_cp = net.shape[2] - 4
        n_row = int(net.shape[0])
        n_col = int(net.shape[1])
        print("network {} shape: ".format(i), net.shape)
        n_points = n_row * n_col
        n_cells = (n_row - 1) * (n_col - 1)
        vtu += "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n".format(n_points, n_cells)

        # format the agps data
        points = str()  # coordinate of each point (x, y, z)
        cells = str()  # ids of each quadrilateral cell (e.g. (0, n_col, n_col + 1, 1) for first cell)
        base_square = np.array((0, n_col, n_col+1, 1), dtype=int)
        for j in range(n_row):
            for k in range(n_col):
                point = net[j,k]
                # add coordinate of a point
                points += "{0} {1} {2}\n".format(point[1], point[2], point[3])
                # add ids of a cell
                if not j == n_row-1 and not k == n_col-1:
                    square =  base_square + (j * n_col + k)
                    square = (str(p) for p in square)
                    cells += " ".join(square) + "\n"

        # add formatted agps data to vtu
        vtu += "      <PointData Scalars=\"scalars\">\n"
        # add point_data to vtu
        for l in range(n_cp):
            vtu += "        <DataArray type=\"Float32\" Name=\"cp{}\" format=\"ascii\">\n".format(l + 1)
            vtu += " ".join(str(cp) for cp in net[:, :, 4 + l].ravel()) + "\n"
            vtu += "        </DataArray>\n"
        vtu += "      </PointData>\n"

        # add points to vtu
        vtu += "      <Points>\n"
        vtu += "        <DataArray type=\"Float32\" Name=\"network{}\" NumberOfComponents=\"3\" " \
               "format=\"ascii\">\n".format(i + 1)
        vtu += points
        vtu += "        </DataArray>\n"
        vtu += "      </Points>\n"

        # add cells to vtu
        vtu += "      <Cells>\n"
        vtu += "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
        vtu += cells
        vtu += "        </DataArray>\n"
        vtu += "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
        vtu += " ".join(str(4*(icell+1)) for icell in range(n_cells)) + "\n"
        vtu += "        </DataArray>\n"
        vtu += "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n"
        vtu += " ".join(str(9) for _ in range(n_cells)) + "\n"
        vtu += "        </DataArray>\n"
        vtu += "      </Cells>\n"

        vtu += "    </Piece>\n"
        vtu += "  </UnstructuredGrid>\n</VTKFile>\n"

        with open(vtu_path, "w") as f:
            f.write(vtu)
    vtm += "  </vtkMultiBlockDataSet>\n</VTKFile>"
    with open("{}.vtm".format(outputname), "w") as f:
        f.write(vtm)


def write_tec(n_wake=0, outputname="agps", inputfile="agps"):
    """convert agps networks to tecplot finite element quadrilaterals"""
    data = read_agps(inputfile) # read agps file & specify the number of networks to omit
    print("n_wake = ", n_wake)
    # write header
    n_headers = data[0].shape[2] # number of headers (e.g. 8 for "irow, x, y, z, cp1, cp2, cp3, cp4")
    n_cp = n_headers - 4 # number of different cps in agps file
    tec = "TITLE = \"AGPS 3D Finite Element Data\"\n"
    tec += "VARIABLES = \"x\", \"y\", \"z\""
    for i in range(n_cp):
        tec += ", \"cp{}\"".format(i+1)
    tec += "\n"

    # write each network as a block
    for i in range(len(data) - n_wake):
        # write the header of the block
        net = data[i]
        n_row = int(net.shape[0])
        n_col = int(net.shape[1])
        print("network {} shape: ".format(i+1), net.shape)
        n_points = n_row * n_col
        n_elements = (n_row-1) * (n_col-1)
        tec += "ZONE T=\"MIXED\", N={}, E={}, DATAPACKING=BLOCK," \
                  " ZONETYPE=FEQUADRILATERAL\n".format(n_points, n_elements)
        # write coordinates (x, y, z) and cps (cp1, cp2, cp3, cp4) in each row
        for l in range(1, n_headers):
            element = net[:, :, l]
            tec += " ".join(map(str, element.ravel())) + "\n"
        # write the ids of each quadrilateral (e.g. (0, n_col, n_col + 1, 1) for first quadrilateral)
        base_square = np.array((0, n_col, n_col + 1, 1)) + 1
        # quads = str()
        # for j in range(n_row-1):
        #     for k in range(n_col-1):
        #         square =  base_square + (j * n_col + k)
        #         square = (str(p) for p in square)
        #         quads += " ".join(square) + "\n"
        # same as the above code, but faster evaluation
        quads = "\n".join("\n".join((" ".join((str(p) for p in (base_square + j * n_col + k))))
                     for k in range(n_col-1))
                    for j in range(n_row-1))
        tec += quads
    with open("{}.dat".format(outputname), "w") as f:
        f.write(tec)