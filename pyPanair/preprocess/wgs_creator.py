#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd


class LaWGS:
    """ A class that defines a geometry in the Langley Wireframe Geometry Standard (LaWGS) format
    The LaWGS geometry is expressed by a list of Networks

    reference: Craidon, C. B., "A Description of the Langley Wireframe Geometry Standard (LaWgs) Format,"
               NASA Technical Memorandum 85767, February 1985.
               Carmichael, R., "LaWGS," Public Domain Aeronautical Software, http://www.pdas.com/lawgs.html,
               [retrieved 19 February 2017]
    """

    def __init__(self, name):
        self.name = name
        self._network_names = list()
        self._networks = list()
        self._boundary_types = list()

    def append_network(self, name, network, boun_type):
        """  append a network to LaWGS
        :param name: name of the network
        :param network: a Network object
        :param boun_type: the boundary type of the Network
        """
        self._network_names.append(name)
        self._networks.append(network)
        self._boundary_types.append(boun_type)

    def create_wgs(self, filename=None):
        """ create a .wgs file from a LaWGS object"""
        if filename is None:
            filename = "{}.wgs".format(self.name)
        wgs = "{} created by wgs_creator\n".format(self.name)
        for (net_name, net) in zip(self._network_names, self._networks):
            wgs += "{}\n".format(net_name)
            num_rows, num_columns = net.shape[:2]
            wgs += "1 {} {} 0   0 0 0   0 0 0    1 1 1  0\n".format(num_rows, num_columns)
            wgs += net.network2wgs()
        with open(filename, "w") as f:
            f.write(wgs)

    def create_aux(self, alpha, mach, cbar, span, sref, xref, zref, filename=None, wgs_filename=None):
        """ create a .aux file (input file for panin) from a LaWGS object"""
        if filename is None:
            filename = "{}.aux".format(self.name)
        if wgs_filename is None:
            wgs_filename = "{}.wgs".format(self.name)
        try:
            alpha = " ".join(map(str, alpha))
        except TypeError:
            if type(alpha) in (int, float):
                pass
            else:
                raise
        boun = " ".join(map(str, self._boundary_types))
        aux = ["// Auxiliary files for {}".format(self.name),
               "WGS {}".format(wgs_filename),
               "MACH {}".format(mach),
               "ALPHA {}".format(alpha),
               "cbar {}".format(cbar),
               "span {}".format(span),
               "sref {}".format(sref),
               "xref {}".format(xref),
               "zref {}".format(zref),
               "BOUN {}".format(boun)]
        aux = "\n".join(aux)
        with open(filename, "w") as f:
            f.write(aux)

    def create_stl(self, filename=None, include_wake=False):
        """ convert LaWGS to stl
        quadrilaterals of the Networks are converted to a pair of trilaterals
        """
        if filename is None:
            filename = "{}.stl".format(self.name)

        def tri2stl(p0, p1, p2):
            """convert coordinates of a trilateral in to a stl facet"""
            norm_vec = np.cross(p1 - p0, p2 - p0)
            stl = "facet normal {0} {1} {2}\n".format(*norm_vec)
            stl += " outer loop\n"
            stl += "  vertex {0} {1} {2}\n".format(*p0)
            stl += "  vertex {0} {1} {2}\n".format(*p1)
            stl += "  vertex {0} {1} {2}\n".format(*p2)
            stl += " endloop\n"
            stl += "endfacet\n"
            return stl

        def quad2stl(p0, p1, p2, p3):
            """convert a quadrilateral into stl facets (or facet if the quadrilateral is actually a trilateral)"""
            if np.array_equal(p0, p1) or np.array_equal(p0, p2):
                return tri2stl(p3, p2, p1)
            elif np.array_equal(p1, p3) or np.array_equal(p2, p3):
                return tri2stl(p0, p1, p2)
            else:
                return tri2stl(p0, p1, p2) + tri2stl(p3, p2, p1)

        stl_all = "solid {}\n".format(self.name)
        for (net, boun) in zip(self._networks, self._boundary_types):
            """ shape of quadrilateral: p0 --- p2
                                        |      |
                                        p1 --- p3"""
            if boun not in (18, 19, 20) or include_wake:
                m = (net.shape[0] - 1) * (net.shape[1] - 1)
                p0 = net[:-1, :-1].reshape((m, 3))
                p1 = net[1:, :-1].reshape((m, 3))
                p2 = net[:-1, 1:].reshape((m, 3))
                p3 = net[1:, 1:].reshape((m, 3))
                stl = (quad2stl(*pnts) for pnts in zip(p0, p1, p2, p3))
                stl = "".join(stl)
                stl_all += stl
        stl_all += "endsolid {}\n".format(self.name)
        with open(filename, "w") as f:
            f.write(stl_all)

    def count_total_panels(self):
        return np.sum(((net.shape[0] - 1) * (net.shape[1] - 1) for net in self._networks))


def read_wgs(filename, wgsname=None, boun_cond=None):
    """ read a LaWGS file and create a LaWGS object
    :param filename: filename of the LaWGS file
    :param wgsname: "name" for the LaWGS object
    :param boun_cond: a tuple containing the boundary conditions of each network
    :return: LaWGS object
    """
    if wgsname is None:
        wgsname = filename.replace(".wgs", "")
    lawgs = LaWGS(wgsname)
    with open(filename, "r") as f:
        f.readline()  # skip first line
        net_id = 0
        while True:
            net_id += 1
            netname = f.readline().split()  # parse the name of the network
            if not netname:
                break
            # parse the header of the network
            header = f.readline().split()
            n_line = int(header[1])
            n_pnt = int(header[2])
            n_total_pnts = 3 * n_line * n_pnt
            cnt = 0
            coords = list()
            while True:
                line = f.readline().split()
                coords += line
                cnt += len(line)
                if cnt == n_total_pnts:  # break after reading all the points in the network
                    break
            network = Network(np.array(coords, dtype=float).reshape((n_line, n_pnt, 3)))
            # set boundary conditions if they are given
            if boun_cond:
                boun = boun_cond[net_id - 1]
            else:
                boun = 1
            lawgs.append_network(netname, network, boun)
        return lawgs


class BasicGeom(np.ndarray):
    _base_shape = (3,)

    def shift(self, distance, inplace=False):
        """shift all the points of a geometry
        :param distance: the distance that points are shifted (e.g. (1., 2., 0.))
        :param inplace: a new object will be created unless inplace=True"""
        distance = np.asarray(distance).reshape(self._base_shape)
        if inplace:
            self[:] += distance
            return self
        else:
            return self + distance

    def _rotate(self, rotmat, rotcenter):
        raise NotImplementedError

    def _rotxyz(self, axis, rotcenter, angle, degrees):
        rotcenter = Point(rotcenter)
        if degrees:
            angle = np.radians(angle)
        if axis == "x":
            rotmat = np.array([[1, 0, 0],
                               [0, np.cos(angle), -np.sin(angle)],
                               [0, np.sin(angle), np.cos(angle)]])
        elif axis == "y":
            rotmat = np.array([[np.cos(angle), 0, np.sin(angle)],
                               [0, 1, 0],
                               [-np.sin(angle), 0, np.cos(angle)]])
        elif axis == "z":
            rotmat = np.array([[np.cos(angle), -np.sin(angle), 0],
                               [np.sin(angle), np.cos(angle), 0],
                               [0, 0, 1]])
        else:
            raise ValueError("axis must be x, y, or z")
        return self._rotate(rotmat, rotcenter)

    def rotx(self, rotcenter, angle, degrees=True):
        return self._rotxyz("x", rotcenter, angle, degrees)

    def roty(self, rotcenter, angle, degrees=True):
        return self._rotxyz("y", rotcenter, angle, degrees)

    def rotz(self, rotcenter, angle, degrees=True):
        return self._rotxyz("z", rotcenter, angle, degrees)


class Network(BasicGeom):
    """ A class that represents a 3 dimensional network which is composed of Lines
    coordinate data is stored in a numpy.ndarray with a shape of (n_lines, n_points, 3)
    """
    _base_shape = (1, 1, 3)

    def __new__(cls, input_network):
        obj = np.asarray(input_network, dtype=float).view(cls)
        if obj.ndim != 3:
            raise ValueError("Network must be 3d array")
        if int(obj.shape[0]) == 1:
            raise ValueError("Network must be composed of at least 2 lines")
        if int(obj.shape[1]) == 1:
            raise ValueError("Line must be composed of at least 2 points")
        if int(obj.shape[2]) != 3:
            raise ValueError("Network must be composed of 3 dimensional points")
        return obj

    def edge(self, edge_number):
        """return the edge of the Network as a Line
        use the plot_wireframe method to check which edge_number corresponds to which edge"""
        edge_number = int(edge_number)
        if edge_number == 1:
            return Line(self[:, 0, :])
        elif edge_number == 2:
            return Line(self[-1])
        elif edge_number == 3:
            return Line(self[:, -1, :])
        elif edge_number == 4:
            return Line(self[0])
        else:
            raise ValueError("edge_number must be 1, 2, 3 or 4")

    def network2wgs(self):
        return "\n".join(("\n".join((" ".join(map("{:.7e}".format, point))
                                     for point in line)) for line in self)) + "\n"

    def concat_row(self, networks):
        """ join edge2 and edge4 of self and input networks
        :param networks: a single Network or a tuple of Networks"""
        if type(networks) is Network:
            networks = (self, networks[1:])
        else:
            if all(type(net) is Network for net in networks):
                networks = [net[1:] for net in networks]
                networks.insert(0, self)
            else:
                raise ValueError("only Networks can be concatenated")
        return Network(np.concatenate(networks))

    def concat_column(self, networks):
        """ join edge1 and edge3 of self and input networks
        :param networks: a single Network or a tuple of Networks"""
        if type(networks) is Network:
            networks = (self, networks[:, 1:])
        else:
            if all(type(net) is Network for net in networks):
                networks = [net[:, 1:] for net in networks]
                networks.insert(0, self)
            else:
                raise ValueError("only Networks can be concatenated")
        return Network(np.concatenate(networks, axis=1))

    def plot_wireframe(self, show_corners=True, show_edges=True, show_normvec=True):
        """ plot the Network as a wireframe
        :param show_corners: display the corner numbers of the Network
        :param show_edges: display the edge numbers of the Network
        :param show_normvec: show a vector pointing out from the front side of the Network
                            (will not work for skewed Networks and Networks consisting only a pair of Lines)"""
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_aspect("equal")
        x, y, z = (self[:, :, i] for i in range(3))
        ax.plot_wireframe(x, y, z)
        # place invisible markers to set the aspect ratio of each axis to be equal
        max_range = np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max()
        xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].ravel() + 0.5 * (x.max() + x.min())
        yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].ravel() + 0.5 * (y.max() + y.min())
        zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].ravel() + 0.5 * (z.max() + z.min())
        for xxb, yyb, zzb in zip(xb, yb, zb):
            ax.plot([xxb], [yyb], [zzb], 'w')
        ax.set_xlabel("$x$", fontsize=16)
        ax.set_ylabel("$y$", fontsize=16)
        ax.set_zlabel("$z$", fontsize=16)
        rowmid, colmid = (i // 2 for i in self.shape[:2])
        if show_corners:
            ax.text(*self[0, 0, :], s="1")
            ax.text(*self[-1, 0, :], s="2")
            ax.text(*self[-1, -1, :], s="3")
            ax.text(*self[0, -1, :], s="4")
        if show_edges:
            ax.text(*self[rowmid, 0, :], s="edge1")
            ax.text(*self[-1, colmid, :], s="edge2")
            ax.text(*self[rowmid, -1, :], s="edge3")
            ax.text(*self[0, colmid, :], s="edge4")
        if show_normvec:
            v_edge4 = self[rowmid, colmid, :] - self[rowmid, colmid - 1, :]
            v_edge1 = self[rowmid, colmid, :] - self[rowmid - 1, colmid, :]
            normvec = np.cross(v_edge1, v_edge4)
            normvec = normvec / np.linalg.norm(normvec) * max_range / 10.
            a = ([pos, vec + pos] for (vec, pos) in zip(normvec, self[rowmid, colmid, :]))
            a = Arrow3D(*a, color="k", mutation_scale=20, arrowstyle="-|>")
            ax.add_artist(a)
        plt.show()

    def trans(self):
        """transpose the row (axis0) and column(axis1) of the Network"""
        return self.transpose((1, 0, 2))

    def _rotate(self, rotmat, rotcenter):
        return Network([np.dot(rotmat, row.T).T for row in self.shift(-rotcenter)]).shift(rotcenter)

    def make_wake(self, edge_number, wake_length):
        """return a wake Network that starts from a edge of the Network
        :param edge_number: the edge number that the wake starts from
        :param wake_length: the length of the wake"""
        edge = self.edge(edge_number)
        wake_end = edge.replace(x=wake_length)
        wake = Network((edge, wake_end))
        wake = wake.trans()
        return wake


class Line(BasicGeom):
    """ A class that represents a 3 dimensional line which is composed of Points
    coordinate data is stored in a numpy.ndarray with a shape of (n_points, 3)
    """
    _base_shape = (1, 3)

    def __new__(cls, input_line):
        obj = np.asarray(input_line, dtype=float).view(cls)
        if obj.ndim != 2:
            raise ValueError("Line must be 2d array")
        if obj.shape[0] == 1:
            raise ValueError("Line must be composed of at least 2 points")
        if obj.shape[1] != 3:
            raise ValueError("Line must be composed of 3 dimensional points")
        return obj

    def linspace(self, stop, num, **kwargs):
        """ create a Network by linearly interpolating the Lines "self" and "stop"
        :param stop: the Line "stop" will become the edge of the interpolated Network
        :param num: number of lines in the interpolated network
        """
        stop = Line(stop)
        if not self.shape == stop.shape:
            raise ValueError("Lines \"self\" and \"stop\" must have the same shape")
        spacing = np.linspace(0., 1., num, **kwargs)
        return Network(self + (stop - self) * spacing[:, np.newaxis, np.newaxis])

    def cosspace(self, stop, num, **kwargs):
        """ create a Network by interpolating the Lines "self" and "stop" using half-cosine spacing
        :param stop: the Line "stop" will become the edge of the interpolated Network
        :param num: number of lines in the interpolated network
        """
        stop = Line(stop)
        if not self.shape == stop.shape:
            raise ValueError("Lines \"self\" and \"stop\" must have the same shape")
        spacing = cosspace(0., 1., num, **kwargs)
        return Network(self + (stop - self) * spacing[:, np.newaxis, np.newaxis])

    def replace(self, x=None, y=None, z=None):
        """ return a Line which has the coordinates of each point replaced by "x", "y", "z"
            coordinates will not be replaced if the variable is None"""
        new_line = [Point(point).replace(x, y, z) for point in self]
        return Line(new_line)

    def flip(self):
        """ return a Line with the order of points flipped from the original"""
        return np.flipud(self)

    def concat(self, lines):
        """concatenate two or more Lines"""
        if type(lines) is Line:
            lines = (self, lines[1:])
        else:
            if all(type(line) is Line for line in lines):
                lines = [line[1:] for line in lines]
                lines.insert(0, self)
            else:
                raise ValueError("only Lines can be concatenated")
        return Line(np.concatenate(lines))

    def _rotate(self, rotmat, rotcenter):
        shift_line = Line(self - rotcenter)
        return Line((rotmat.dot(shift_line.T)).T + rotcenter)

    def split_half(self):
        """ split the Line into two halves and return copies of the halves
        useful when splitting an airfoil into upper and lower Lines"""
        first_half = Line(np.array(self[:self.shape[0] // 2 + 1]))
        second_half = Line(np.array(self[self.shape[0] // 2:]))
        return first_half, second_half


def read_airfoil(filename, span_pos=0., expansion_ratio=1.):
    """ read the coordinates of a airfoil from a csv file and create Line from it
    see the examples "naca2412.csv" and "rae2822.csv" in the examples directory
    upper and lower surfaces must have the same number of points
    first and last points of the upper and lower surfaces must coincide
    :param filename: name of the csv file
    :param span_pos: the y-axis coordinate will be set at span_pos
    :param expansion_ratio: the chord length and thickness (i.e. x,z-axis coordinates) will be multiplied by this ratio
    """
    afoil = pd.read_csv(filename)
    # assert that the upper and lower surfaces have the same number of points
    if np.unique(afoil.count()).size >= 2:
        raise ValueError("upper and lower surfaces must have the same number of points")
    # assert that the first and last points of the upper and lower surfaces coincide
    if not np.array_equal(afoil.head(1)[["xup", "zup"]], afoil.head(1)[["xlow", "zlow"]]):
        raise ValueError("first points of the upper and lower surfaces must coincide")
    if not np.array_equal(afoil.tail(1)[["xup", "zup"]], afoil.tail(1)[["xlow", "zlow"]]):
        raise ValueError("last points of the upper and lower surfaces must coincide")
    # convert the airfoil csv into a Line
    xup = np.flipud(afoil["xup"])
    zup = np.flipud(afoil["zup"])
    xlow = afoil["xlow"].values
    zlow = afoil["zlow"].values
    n_pnts = xup.shape[0]
    afoil_line = np.ones((n_pnts * 2 - 1, 3))
    afoil_line[:, 1] *= span_pos
    afoil_line[:n_pnts, 0] = xup
    afoil_line[:n_pnts, 2] = zup
    afoil_line[n_pnts:, 0] = xlow[1:]
    afoil_line[n_pnts:, 2] = zlow[1:]
    afoil_line[:, [0, 2]] *= expansion_ratio
    return Line(afoil_line)


def naca4digit(digits, num=25, chord=1., span_pos=0.):
    """ create a Line of a naca 4 digit airfoil
    :param digits: the 4 digits that specify the airfoil (e.g. 2412)
    :param num: number of points for upper and lower surfaces
    :param chord: chord length of airfoil
    :param span_pos: y-coordinate of the Line
    """
    digits = str(digits)
    if len(digits) != 4:
        raise ValueError("enter 4 digits (e.g. 2412)")
    m = float(digits[0]) / 100.
    p = float(digits[1]) / 10.
    t = float(digits[2:]) / 100.0
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # originally -0.1015, but -0.1036 for sharp edge
    x = cosspace(0, 1, num)
    zt = 5 * t * (a0 * (x ** 0.5) + a1 * x + a2 * (x ** 2) + a3 * (x ** 3) + a4 * (x ** 4))  # thickness

    def camber_gradient(xi, p):
        if 0 <= xi < p:
            zc = (m / (p ** 2)) * (2 * p * xi - xi ** 2)  # camber
            grad_zc = (2 * m / (p ** 2)) * (p - xi)  # gradient
        else:
            zc = (m / ((1 - p) ** 2)) * (1 - 2 * p + 2 * p * xi - xi ** 2)  # camber
            grad_zc = (2 * m / ((1 - p) ** 2)) * (p - xi)  # gradient
        return zc, grad_zc

    vfunc = np.vectorize(camber_gradient)
    zc, grad_zc = vfunc(x, p)
    theta = np.arctan(grad_zc)
    xup = chord * (x - zt * np.sin(theta))
    zup = chord * (zc + zt * np.cos(theta))
    xup = np.flipud(xup)
    zup = np.flipud(zup)
    xlow = chord * (x + zt * np.sin(theta))
    zlow = chord * (zc - zt * np.cos(theta))
    # modify the zup & zlow to make the thickness 0 at the TE
    zup[0] = 0.
    zlow[-1] = 0.
    afoil_line = np.ones((num * 2 - 1, 3))
    afoil_line *= span_pos
    afoil_line[:num, 0] = xup
    afoil_line[:num, 2] = zup
    afoil_line[num:, 0] = xlow[1:]
    afoil_line[num:, 2] = zlow[1:]
    return Line(afoil_line)


class Point(BasicGeom):
    """ A class that represents a 3 dimensional point
    coordinate data is stored in a numpy.ndarray with a shape of (3,)
    """
    _base_shape = (3,)

    def __new__(cls, input_point):
        obj = np.asarray(input_point, dtype=float).view(cls)
        if obj.size != 3:
            raise ValueError("Point must be 3 dimensional")
        return obj.reshape(cls._base_shape)

    def linspace(self, stop, num, **kwargs):
        """ create a Line by linearly interpolating the Points "self" and "stop"
        :param stop: the Point "stop" will become the end of the interpolated Line
        :param num: number of points in the interpolated line
        """
        stop = Point(stop)
        spacing = np.linspace(0., 1., num, **kwargs)
        return Line(self + (stop - self) * spacing[:, np.newaxis])

    def cosspace(self, stop, num, **kwargs):
        """ create a Line by interpolating the Points "self" and "stop" using half-cosine spacing
        :param stop: the Point "stop" will become the end of the interpolated Line
        :param num: number of points in the interpolated line
        """
        stop = Point(stop)
        spacing = cosspace(0., 1., num, **kwargs)
        return Line(self + (stop - self) * spacing[:, np.newaxis])

    def replace(self, x=None, y=None, z=None):
        """ return a point which has coordinates replaced by "x", "y", "z"
            coordinates will not be replaced if the variable is None"""
        new_point = [s if c is None else c for (s, c) in zip(self, (x, y, z))]
        return Point(new_point)

    def square_network(self, p2, p3, p4, row_num, column_num, row_spacing="linspace", column_spacing="linspace",
                       **kwargs):
        """ create a Network by interpolating the Points self, p2, p3, and p4
        the columns will become edge 1 and 3
        the rows will become edge 2 and 4
         p4  edge3  p3

        edge4      edge2

         p1  edge1  p2"""
        if row_spacing in ("linspace", "lin"):
            line1 = self.linspace(p4, row_num, **kwargs)
            line2 = p2.linspace(p3, row_num, **kwargs)
        elif row_spacing in ("cosspace", "cos"):
            line1 = self.cosspace(p4, row_num, **kwargs)
            line2 = p2.cosspace(p3, row_num, **kwargs)
        else:
            raise ValueError("row_spacing must be linspace or cosspace")
        if column_spacing in ("linspace", "lin"):
            return line1.linspace(line2, column_num, **kwargs)
        elif column_spacing in ("cosspace", "cos"):
            return line1.cosspace(line2, column_num, **kwargs)
        else:
            raise ValueError("spacing must be linspace or cosspace")


def cosspace(start, stop, num, **kwargs):
    """ half cosine interpolation between "start" and "stop"
    """
    return start + (stop - start) * 0.5 * (1 - np.cos(np.linspace(0., np.pi, num, **kwargs)))


class Arrow3D(FancyArrowPatch):
    """a class for drawing a 3D arrow"""

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)
