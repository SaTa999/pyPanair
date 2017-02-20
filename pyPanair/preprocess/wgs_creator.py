#!/usr/bin/env python
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from pyPanair.utilities import cosspace, naca4digit, bspline


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
        wgs = "{} created from wgs_creator\n".format(self.name)
        for (net_name, net) in zip(self._network_names, self._networks):
            wgs += "{}\n".format(net_name)
            num_rows, num_columns = net.shape[:2]
            wgs += "1 {} {} 0   0 0 0   0 0 0    1 1 1  0\n".format(num_rows, num_columns)
            wgs += net.network2wgs()
        with open(filename, "w") as f:
            f.write(wgs)

    def create_aux(self, filename, alpha, mach, cbar, span, sref, xref, zref, wgs_filename=None):
        """ create a .aux file (input file for panin) from a LaWGS object"""
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
        aux = ["// PAN AIR CASE MANUAL CASE 1",
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

    def create_stl(self, filename, include_wake=False):
        """ convert LaWGS to stl
        quadrilaterals of the Networks are converted to a pair of trilaterals
        """
        def tri2stl(p0, p1, p2):
            """convert coordinates of a trilateral in to a stl facet"""
            norm_vec = np.cross(p1-p0, p2-p0)
            stl =  "facet normal {0} {1} {2}\n".format(*norm_vec)
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
            if not boun in (18, 19, 20) or include_wake:
                m = (net.shape[0] - 1) * (net.shape[1] - 1)
                p0 = net[:-1,:-1].reshape((m, 3))
                p1 = net[1:,:-1].reshape((m, 3))
                p2 = net[:-1,1:].reshape((m, 3))
                p3 = net[1:,1:].reshape((m, 3))
                stl = (quad2stl(*pnts)  for pnts in zip(p0,p1,p2,p3))
                stl = "".join(stl)
                stl_all += stl
        stl_all += "endsolid {}\n".format(self.name)
        with open(filename, "w") as f:
            f.write(stl_all)

    def count_total_panels(self):
        return np.sum(((net.shape[0]-1)*(net.shape[1]-1) for net in self._networks))


def read_wgs(filename, wgsname=None, boun_cond=None):
    """ read a LaWGS file and create a LaWGS object
    :param filename: filename of the LaWGS file
    :param wgsname: "name" for the LaWGS object
    :param bonu_cond: a tuple containing the boundary conditions of each network
    :return: LaWGS object
    """
    if wgsname is None:
        wgsname = filename.replace(".wgs", "")
    lawgs = LaWGS(wgsname)
    with open(filename, "r") as f:
        f.readline() # skip first line
        net_id = 0
        while True:
            net_id += 1
            netname = f.readline().split() # parse the name of the network
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
                if cnt == n_total_pnts: # break after reading all the points in the network
                    break
            network = np.array(coords, dtype=float).reshape((n_line, n_pnt, 3))
            # set boundary conditions if they are given
            if boun_cond:
                boun = boun_cond[net_id]
            else:
                boun = 1
            lawgs.append_network(netname, network, boun)
        return lawgs


class Network(np.ndarray):
    """
    wgs形式で使用するネットワークのクラス
    Lineの集合として表す
    形式はndim=3のnumpy-ndarray
    """

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
        """NetworkのedgeをLineとして返す"""
        edge_number = int(edge_number)
        if edge_number is 1:
            return Line(self[:, 0, :])
        elif edge_number is 2:
            return Line(self[-1])
        elif edge_number is 3:
            return Line(self[:, -1, :])
        elif edge_number is 4:
            return Line(self[0])
        else:
            raise ValueError("edge_number must be 1, 2, 3 or 4")

    def network2wgs(self):
        return "\n".join(("\n".join((" ".join(map("{:.7e}".format, point))
                          for point in line)) for line in self)) + "\n"

    def concat_row(self, networks):
        """Networkのedge2,4同士を繋げたNetworkを返す"""
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
        """Networkのedge1,3同士を繋げたNetworkを返す"""
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
        """Networkをワイヤーフレームとして表示する"""
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_aspect("equal")
        X, Y, Z = (self[:, :, i] for i in range(3))
        ax.plot_wireframe(X, Y, Z)
        # 各方向のアスペクト比を同じにするため隅に透明なマーカーを置く
        max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max()
        Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].ravel() + 0.5 * (X.max() + X.min())
        Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].ravel() + 0.5 * (Y.max() + Y.min())
        Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].ravel() + 0.5 * (Z.max() + Z.min())
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')
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

    def shift(self, shift):
        return self + shift

    def trans(self):
        """Networkの行と列を入れ替える"""
        return self.transpose((1, 0, 2))

    def rotz(self, center, degrees):
        center = Point(center)
        rad = np.radians(degrees)
        rz = np.array([[np.cos(rad), -np.sin(rad), 0],
                       [np.sin(rad), np.cos(rad), 0],
                       [0, 0, 1]])
        return Network([np.dot(rz, row.T).T for row in self.shift(-center)]).shift(center)

    def make_wake(self, edge_number, wake_length):
        """Networkのedgeを始点とするwakeのNetworkを返す"""
        edge = self.edge(edge_number)
        wake_end = edge.shifted_line(x=wake_length)
        wake = Network((edge, wake_end))
        wake = wake.trans()
        return wake


class Line(np.ndarray):
    """
    Networkを構成するLineのクラス
    Pointの集合として表す
    形式はndim=2のnumpy-ndarray
    """

    def __new__(cls, input_line):
        obj = np.asarray(input_line, dtype=float).view(cls)
        if obj.ndim is not 2:
            raise ValueError("Line must be 2d array")
        if int(obj.shape[0]) is 1:
            raise ValueError("Line must be composed of at least 2 points")
        if int(obj.shape[1]) is not 3:
            raise ValueError("Line must be composed of 3 dimensional points")
        return obj

    def linspace(self, stop, num, **kwargs):
        """自身からstopまでをnum分割したLineを返す（等間隔）"""
        if not self.shape[0] == stop.shape[0]:
            raise ValueError("length of lines are unequal")
        lin = np.linspace(0., 1., num, **kwargs)
        return Network(self + (stop - self) * lin[:, np.newaxis, np.newaxis])

    def cosspace(self, stop, num, **kwargs):
        """自身からstopまでをnum分割したLineを返す（half-cosine spacing）"""
        if not self.shape[0] == stop.shape[0]:
            raise ValueError("length of lines are unequal")
        cs = cosspace(0., 1., num, **kwargs)
        return Network(self + (stop - self) * cs[:, np.newaxis, np.newaxis])

    def shifted_line(self, x=None, y=None, z=None):
        """Lineの全Pointを指定の座標を変更した新しいLineを返す"""
        new_line = [Point(point).shifted_point(x, y, z) for point in self]
        return Line(new_line)

    def flip(self):
        """向きを反転したLineを返す"""
        return np.flipud(self)

    def concat(self, lines):
        """Lineを繋げたLineを返す"""
        if type(lines) is Line:
            lines = (self, lines[1:])
        else:
            if all(type(line) is Line for line in lines):
                lines = [line[1:] for line in lines]
                lines.insert(0, self)
            else:
                raise ValueError("only Lines can be concatenated")
        return Line(np.concatenate(lines))

    def rotx(self, center, degrees):
        center = Point(center)
        rad = np.radians(degrees)
        rz = np.array([[1, 0, 0],
                       [0, np.cos(rad), -np.sin(rad)],
                       [0, np.sin(rad), np.cos(rad)]])
        shift_line = Line(self - center)
        return Line((rz @ shift_line.T).T + center)


class Point(np.ndarray):
    """
    Lineを構成するPointのクラス
    中身はnumpyのndarray
    """

    def __new__(cls, input_point):
        obj = np.asarray(input_point, dtype=float).view(cls)
        if obj.ndim is not 1:
            raise ValueError("Point must be 1d array")
        if int(obj.shape[0]) is not 3:
            raise ValueError("Point must be 3 dimensional")
        return obj

    def linspace(self, stop, num, **kwargs):
        """自身からstopまでをnum分割したLineを返す（等間隔）"""
        lin = np.linspace(0., 1., num, **kwargs)
        return Line(self + (stop - self) * lin[:, np.newaxis])

    def cosspace(self, stop, num, **kwargs):
        """自身からstopまでをnum分割したLineを返す（half-cosine spacing）"""
        cs = cosspace(0., 1., num, **kwargs)
        return Line(self + (stop - self) * cs[:, np.newaxis])

    def shifted_point(self, x=None, y=None, z=None):
        """Pointの指定の座標を変更した新しいPointを返す"""
        new_point = [self[i] if c is None else c for (i, c) in enumerate((x, y, z))]
        return Point(new_point)

    def square_network(self, p2, p3, p4, row_num, column_num, row_spacing="linspace", column_spacing="linspace",
                       **kwargs):
        """自身及びp2～p4をnetwork cornerとするNetworkを生成する
        edge1,3がcolumnに、edge2,4がrowに対応
        edge2,4：x軸方向（気流方向）、edge1,3：y軸方向とすれば間違いはないはず
        z軸は下図で紙面上向き
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
            return line1.linspace(line2, column_num)
        elif column_spacing in ("cosspace", "cos"):
            return line1.cosspace(line2, column_num)
        else:
            raise ValueError("spacing must be linspace or cosspace")


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def main(xfin=20., yfin=192., zfin=0.5, delta=0., span_fin=50., rchord_fin=80., taper=1., sweep=0.):
    include_fins = True  # 垂直フィンを含める場合True（Falseの場合はクリーン形態）

    wgs_aux = LaWGS("model3")
    # 主翼のwgsを記述する
    # <editor-fold desc="main_wing">
    np_wing_x = 35
    np_wing_y = 35
    np_wing_z = 4

    chord_wing = 80.
    halfspan_wing = 202.
    naca2410u = partial(naca4digit, 2, 4, 10, surface = 'upper')
    naca2410l = partial(naca4digit, 2, 4, 10, surface = 'lower')
    X = cosspace(chord_wing, 0.0, np_wing_x)
    shift = np.array([0., 22., 3.85])

    gridu = [naca2410u(x, chord_wing) for x in X]
    gridl = [naca2410l(x, chord_wing) for x in X]
    gridl.reverse()
    wingrootu = Line([[gr[0], 0., gr[1]] for gr in gridu])
    wingrootu = wingrootu + shift
    wingtipu = wingrootu.shifted_line(y=halfspan_wing)
    wingu = wingrootu.linspace(wingtipu, np_wing_y)
    wingrootl = Line([[gr[0], 0., gr[1]] for gr in gridl])
    wingrootl = wingrootl + shift
    wingtipl = wingrootl.shifted_line(y=halfspan_wing)
    wingl = wingrootl.linspace(wingtipl, np_wing_y)
    wing = wingu.concat_column(wingl)
    wgs_aux.append_network("wing", wing, 1)

    wingtip = wingtipu.linspace(wingtipl.flip(), np_wing_z)
    wgs_aux.append_network("wingtip", wingtip, 1)
    # </editor-fold>

    # 舵翼のwgsを記述する
    # <editor-fold desc="fins">
    if include_fins:
        # rchord_fin = 80.
        # span_fin = 50.
        # taper = 1.
        # delta = 0
        # xfin = 0.
        # yfin = 130.
        tchord_fin = rchord_fin * taper

        np_fin_x = 30#int(rchord_fin // 3)
        np_fin_y = 4
        np_fin_z = 25#int(span_fin // 3)

        shift = Point([xfin-0.25*rchord_fin, yfin, 3.85 + 5.54 + zfin])
        center = Point([rchord_fin * 0.25, 0., 0.])
        rtshift = Point([(rchord_fin - tchord_fin) * 0.25+span_fin*np.sin(np.radians(sweep)), 0., 0.])
        finfoilu = partial(naca4digit, 0, 0, 10, surface='upper')
        finfoill = partial(naca4digit, 0, 0, 10, surface='lower')
        rootX = cosspace(rchord_fin, 0.0, np_fin_x)
        tipX = cosspace(tchord_fin, 0.0, np_fin_x)

        grid1 = [finfoill(x, rchord_fin) for x in rootX]
        grid2 = [finfoill(x, tchord_fin) for x in tipX]
        fini1 = Line([[gr[0], gr[1], 0.] for gr in grid1])
        fini2 = Line([[gr[0], gr[1], span_fin] for gr in grid2]) + rtshift
        fini = fini1.linspace(fini2, np_fin_z)
        fini = fini.rotz(center, delta)
        fini = fini.shift(shift)
        wgs_aux.append_network("fininner", fini, 1)

        grid1 = [finfoilu(x, rchord_fin) for x in rootX]
        grid2 = [finfoilu(x, tchord_fin) for x in tipX]
        fino1 = Line([[gr[0], gr[1], 0.] for gr in grid1])
        fino2 = Line([[gr[0], gr[1], span_fin] for gr in grid2]) + rtshift
        fino = fino1.linspace(fino2, np_fin_z)
        fino = Network(np.fliplr(fino))
        fino = fino.rotz(center, delta)
        fino = fino.shift(shift)
        wgs_aux.append_network("finouter", fino, 1)

        fintop1 = fini.edge(2)
        fintop2 = fino.edge(2)
        fintop1 = fintop1.flip()
        fintop = fintop2.linspace(fintop1, np_fin_y)
        wgs_aux.append_network("fintiptop", fintop, 1)

        finbot1 = fini.edge(4)
        finbot2 = fino.edge(4)
        finbot1 = finbot1.flip()
        finbot = finbot1.linspace(finbot2, np_fin_y)
        wgs_aux.append_network("fintipbottom", finbot, 1)
    # </editor-fold>

    # 胴体中央部のwgsを記述する
    # <editor-fold desc="body">
    np_body_y = 5
    np_body_z1 = 5
    np_body_z2 = 3

    body_z1 = 0.
    body_z2 = 55.
    body_y1 = 0.

    body1 = wingu.edge(4).flip()
    body2 = body1.shifted_line(z=body_z2)
    body3 = body2.shifted_line(y=body_y1)
    bodysideupper = body1.cosspace(body2, np_body_z1)
    bodytop = body2.cosspace(body3, np_body_y)
    bodyupper = bodysideupper.concat_row(bodytop)
    wgs_aux.append_network("bodyupper", bodyupper, 1)

    body6 = wingl.edge(4)
    body5 = body6.shifted_line(z=body_z1)
    body4 = body5.shifted_line(y=body_y1)
    bodybot = body4.cosspace(body5, np_body_y)
    bodysidelower = body5.cosspace(body6, np_body_z2)
    bodylower = bodybot.concat_row(bodysidelower)
    wgs_aux.append_network("bodylower", bodylower, 1)
    # </editor-fold>

    # noseのwgsを記述する
    # <editor-fold desc="nose">
    np_nose_x1 = 15
    np_nose_x2 = 5
    np_nose_y = 5
    np_nose_z1 = 3
    np_nose_z2 = 5

    nose_x1 = -77.
    nose_x2 = 0.
    nose_y1 = 0.
    nose_y2 = nose_y1 + 12.5
    nose_y3 = 22.
    nose_z1 = 0.
    nose_z2 = nose_z1 + 3.85
    nose_z3 = nose_z1 + 13.
    nose_z5 = nose_z3 + 25.
    nose_z6 = 55.
    nose_z4 = nose_z3 + (nose_z2 - nose_z1) / (nose_z6 - nose_z1) * (nose_z5 - nose_z3)

    cp2 = np.array([[-77.000, 0.000, 25.500],
                    [-77.000, 7.894, 17.606],
                    [-67.570, 13.663, 11.408],
                    [-57.000, 14.968, 9.623]])
    cp4 = np.array([[-77.000, 0.000, 25.500],
                    [-77.000, 8.009, 33.509],
                    [-67.617, 13.658, 40.072],
                    [-57.000, 14.968, 42.416]])
    spl2 = bspline(cp2)
    spl4 = bspline(cp4)
    tcks = cosspace(0, 1, np_nose_x1)
    nose21 = Line(spl2(tcks))
    nose41 = Line(spl4(tcks))
    p2 = Point((nose_x2, nose_y3, nose_z1))
    p4 = p2.shifted_point(z=nose_z6)
    nose22 = Point(nose21[-1]).cosspace(p2, np_nose_x2)
    nose42 = Point(nose41[-1]).cosspace(p4, np_nose_x2)
    nose2 = nose21.concat(nose22)
    nose4 = nose41.concat(nose42)

    nose1 = nose2.shifted_line(y=nose_y1)
    nose5 = nose4.shifted_line(y=nose_y1)
    nose3 = Line((nose2 * (nose_z6 - nose_z2) + nose4 * nose_z2) / nose_z6)

    nosebot = nose1.cosspace(nose2, np_nose_y)
    noselow = nose2.cosspace(nose3, np_nose_z1)
    noseup = nose3.cosspace(nose4, np_nose_z2)
    nosetop = nose4.cosspace(nose5, np_nose_y)
    nose = nosebot.concat_row((noselow, noseup, nosetop))
    wgs_aux.append_network("nose", nose, 1)
    # </editor-fold>

    # 胴体後方部のwgsを記述する
    # <editor-fold desc="tail">
    np_tail_x1 = 4
    np_tail_x2 = 9
    np_tail_y = 5
    np_tail_z1 = 5
    np_tail_z2 = 3
    tail_x1 = 0.
    tail_x2 = tail_x1 + 6.
    tail_x3 = tail_x2 + 216.
    tail_y1 = 0.
    tail_y2 = tail_y1 + 3.
    tail_y4 = tail_y1 + 22.
    tail_y3 = tail_y4 - 0.1
    tail_z1 = 0.
    tail_z2 = 3.85
    tail_z3 = tail_z1 + 40.
    tail_z5 = tail_z1 + 55.
    tail_z4 = tail_z3 + (tail_z2 - tail_z1) / (tail_z5 - tail_z1) * (tail_z5 - tail_z3)
    shift = np.array([80., 0., 0.])

    p1 = Point([tail_x1, tail_y1, tail_z1])
    p2 = p1.shifted_point(y=tail_y4)
    p3 = p2.shifted_point(z=tail_z2)
    p4 = p3.shifted_point(z=tail_z5)
    p5 = p4.shifted_point(y=tail_y1)
    p6 = p5.shifted_point(x=tail_x2)
    p7 = p6.shifted_point(y=tail_y3)
    p8 = p7.shifted_point(z=tail_z2)
    p9 = p8.shifted_point(z=tail_z1)
    p10 = p9.shifted_point(y=tail_y1)
    tail1bot = p1.square_network(p2, p9, p10, np_tail_x1, np_tail_y, "cos", "cos")
    tail1low = p2.square_network(p3, p8, p9, np_tail_x1, np_tail_z1, "cos", "cos")
    tail1up = p3.square_network(p4, p7, p8, np_tail_x1, np_tail_z2, "cos", "cos")
    tail1top = p4.square_network(p5, p6, p7, np_tail_x1, np_tail_y, "cos", "cos")

    tail1lower = tail1bot.concat_row(tail1low)
    tail1upper = tail1up.concat_row(tail1top)

    p1 = p10
    p2 = p9
    p3 = p8
    p4 = p7
    p5 = p6
    p6 = Point([tail_x3, tail_y1, tail_z5])
    p7 = p6.shifted_point(y=tail_y2)
    p8 = p7.shifted_point(z=tail_z4)
    p9 = p8.shifted_point(z=tail_z3)
    p10 = p9.shifted_point(y=tail_y1)
    tail2bot = p1.square_network(p2, p9, p10, np_tail_x1, np_tail_y, "cos", "cos")
    tail2low = p2.square_network(p3, p8, p9, np_tail_x1, np_tail_z1, "cos", "cos")
    tail2up = p3.square_network(p4, p7, p8, np_tail_x1, np_tail_z2, "cos", "cos")
    tail2top = p4.square_network(p5, p6, p7, np_tail_x1, np_tail_y, "cos", "cos")

    tail2lower = tail2bot.concat_row(tail2low)
    tail2upper = tail2up.concat_row(tail2top)

    tailupper = tail1upper.concat_column(tail2upper)
    tailupper = tailupper.shift(shift)
    taillower = tail1lower.concat_column(tail2lower)
    taillower = taillower.shift(shift)

    wgs_aux.append_network("tailupper", tailupper, 1)
    wgs_aux.append_network("taillower", taillower, 1)

    tailbaseu1 = tail2up.edge(3).flip()
    tailbaseu2 = tailbaseu1.shifted_line(y=tail_y1)
    tailbaseu = tailbaseu1.cosspace(tailbaseu2, np_tail_y)
    tailbaseu = tailbaseu.shift(shift)
    wgs_aux.append_network("tailbaseupper", tailbaseu, 5)

    tailbasel1 = tail2low.edge(3).flip()
    tailbasel2 = tailbasel1.shifted_line(y=tail_y1)
    tailbasel = tailbasel1.cosspace(tailbasel2, np_tail_y)
    tailbasel = tailbasel.shift(shift)
    wgs_aux.append_network("tailbaselower", tailbasel, 5)
    # </editor-fold>

    # wakeのwgsを記述する
    # <editor-fold desc="wake">
    wake_length = chord_wing * 25.

    wingwake = wing.make_wake(3, wake_length)
    wgs_aux.append_network("wingwake", wingwake, 18)

    if include_fins:
        finwake = fini.make_wake(1, wake_length)
        wgs_aux.append_network("finwake", finwake, 18)

    bodybaseuwake = tailupper.make_wake(3, wake_length)
    wgs_aux.append_network("bodybaseuwake", bodybaseuwake, 18)

    bodybaselwake = taillower.make_wake(3, wake_length)
    wgs_aux.append_network("bodybaselwake", bodybaselwake, 18)

    bodywake = tailupper.make_wake(4, wake_length)
    wgs_aux.append_network("bodywake", bodywake, 20)
    # </editor-fold>

    wgs_aux.create_wgs()

    #alpha = np.arange(-6, 18, 2)
    #alpha = list(zip(*[iter(alpha)] * 3))
    alpha = (2,)
    for (i, a) in enumerate(alpha):
        wgs_aux.create_aux("model3_{}.aux".format(i + 1), alpha=a, mach=0.0587, cbar=80, span=404, sref=32320,
                           xref=33, zref=30)

    wgs_aux.create_stl("model3.stl")

    delta1 = read_wgs("delta1.wgs")
    delta1.create_stl("delta1.stl")


if __name__ == '__main__':
    main()
