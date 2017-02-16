#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 01 21:52:31 2016

@author: Satoshi
LaWGS形式ファイルを作成するプログラム"""
from math import atan, cos, pi, sin, radians
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import InterpolatedUnivariateSpline, splev


class WGS_AUX(object):
    """wgsとauxを生成するクラス"""

    def __init__(self, name):
        self._name = name
        self._network_names = []
        self._networks = []
        self._boun = []

    def add_network(self, name, network, boun):
        """
        networkを追加する
        （ついでに境界条件（BOUN）を保存しておく）
        """
        self._network_names.append(name)
        self._networks.append(network)
        self._boun.append(boun)

    def create_wgs(self, filename):
        """
        wgsにnetworkを追加する
        （ついでに境界条件（BOUN）を保存しておく）
        """
        wgs = self._name + "\n"
        for (net_name, net) in zip(self._network_names, self._networks):
            wgs += "{0}{1}".format(net_name, "\n")
            num_rows, num_columns = net.shape[:2]
            wgs += "1 {0} {1} 0   0 0 0   0 0 0    1 1 1  0{2}".format(num_rows, num_columns, "\n")
            wgs += net.network2wgs()
        with open(filename, "w") as f:
            f.write(wgs)

    def create_aux(self, filename, aoa):
        """auxを生成する．適宜内容を書き換えること．"""
        try:
            aoa = " ".join(map(str, aoa))
        except TypeError as e:
            if type(aoa) in (int, float):
                pass
            else:
                raise
        boun = " ".join(map(str, self._boun))
        aux = ["// PAN AIR CASE MANUAL CASE 1",
               "WGS model3.wgs",
               "MACH 0.0587",
               "ALPHA {}".format(aoa),
               "sref 32320.0",
               "span 404.0",
               "cbar 80.0",
               "xref 33.0",
               "zref 30.0",
               "BOUN {}".format(boun)]
        aux = "\n".join(aux)
        with open(filename, "w") as f:
            f.write(aux)

    def create_stl(self, filename, include_wake=False):
        """Networkをstl形式に変換、保存する"""
        stl = "solid {}".format(filename)
        stl_all = ["solid model3.stl\n"]
        def net2stl(p0, p1, p2, nv1):
            stl =  "facet normal {0} {1} {2}\n".format(*nv1)
            stl += " outer loop\n"
            stl += "  vertex {0} {1} {2}\n".format(*p0)
            stl += "  vertex {0} {1} {2}\n".format(*p1)
            stl += "  vertex {0} {1} {2}\n".format(*p2)
            stl += " endloop\n"
            stl += "endfacet\n"
            return stl
        
        num_pan = 0
        for (net, boun) in zip(self._networks, self._boun):
            if not boun in (18, 19, 20) or include_wake:
                m = (net.shape[0] - 1) * (net.shape[1] - 1)
                num_pan += m
                p0 = net[:-1,:-1].reshape((m, 3))
                p1 = net[1:,:-1].reshape((m, 3))
                p2 = net[:-1,1:].reshape((m, 3))
                p3 = net[1:,1:].reshape((m, 3))
                nv1 = np.array([np.cross(u,v) for (u,v) in zip(p1-p0,p2-p0)])
                nv2 = np.array([np.cross(u,v) for (u,v) in zip(p2-p3,p1-p3)])
                stl = [net2stl(po0,po1,po2,nov1)+net2stl(po3,po2,po1,nov2)  for (po0, po1, po2, po3, nov1, nov2) in zip(p0,p1,p2,p3,nv1,nv2)]
                stl = "".join(stl)
                stl_all.append(stl)
        print(num_pan)
        stl_all.append("endsolid model3.stl")
        stl_all = "".join(stl_all)
        with open(filename, "w") as f:
            f.write(stl_all)


class Network(np.ndarray):
    """
    wgs形式で使用するネットワークのクラス
    Lineの集合として表す
    形式はndim=3のnumpy-ndarray
    """

    def __new__(cls, input_network):
        obj = np.asarray(input_network, dtype=float).view(cls)
        if obj.ndim is not 3:
            raise ValueError("Network must be 3d array")
        if int(obj.shape[0]) is 1:
            raise ValueError("Network must be composed of at least 2 lines")
        if int(obj.shape[1]) is 1:
            raise ValueError("Line must be composed of at least 2 points")
        if int(obj.shape[2]) is not 3:
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
        return " ".join(["\n".join([" ".join(map(str, point)) for point in line]) for line in self]) + "\n"

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
        rad = radians(degrees)
        rz = np.array([[cos(rad), -sin(rad), 0],
                       [sin(rad), cos(rad), 0],
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


def naca4digit(m, p, xx, x, c, surface):
    """
    naca4digit翼型の2次元座標を返す
    """
    m /= 100.0
    p /= 10.0
    t = xx / 100.0
    xi = x / c
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # 通常は-0.1015 sharp edgeの場合は-0.1036
    yt = 5 * t * (a0 * (xi ** 0.5) + a1 * xi + a2 * (xi ** 2) + a3 * (xi ** 3) + a4 * (xi ** 4))
    if x < 0. or c < x:
        raise ValueError("x = {0} is out of range for c = {1}".format(x, c))
    else:
        pass
    # TEの厚みが0になるように修正
    if x == c:
        return [x, 0.]
    else:
        pass
    if 0 <= xi < p:
        yc = (m / (p ** 2)) * (2 * p * xi - xi ** 2)  # camber
        grad_yc = (2 * m / (p ** 2)) * (p - xi)  # gradient
    else:
        yc = (m / ((1 - p) ** 2)) * (1 - 2 * p + 2 * p * xi - xi ** 2)  # camber
        grad_yc = (2 * m / ((1 - p) ** 2)) * (p - xi)  # gradient
    theta = atan(grad_yc)
    if surface == 'upper':
        xu = c * (xi - yt * sin(theta))
        yu = c * (yc + yt * cos(theta))
        return [xu, yu]
    elif surface == 'lower':
        xu = c * (xi + yt * sin(theta))
        yu = c * (yc - yt * cos(theta))
        return [xu, yu]
    else:
        raise ValueError('Error! surface must be \'upper\' or \'lower\'')


def flatplate(p, t, x, c):
    """
    先端が楕円状の薄板翼型の座標を返す
    """
    b = p * c
    if 0 <= x < b:
        xi = b - x
        y = t * (1. - (xi / b) ** 2) ** 0.5
        return [x, y]
    elif c - b < x <= c:
        xi = b - c + x
        y = t * (1. - (xi / b) ** 2) ** 0.5
        return [x, y]
    elif b <= x <= c - b:
        return [x, t]
    else:
        raise ValueError("x = {0}  out of bounds for c = {1}".format(x, c))


def cosspace(start, stop, num, **kwargs):
    return start + (stop - start) * 0.5 * (1 - np.cos(np.linspace(0., pi, num, **kwargs)))


def bezierspline(cv, degree=3, periodic=False):
    """ ベジエ曲線の関数を返す
        cv :      制御点のarray
        degree:   曲線の次数
        periodic: True - 閉曲線
                  False - 開曲線
    """
    # If periodic, extend the point array by count+degree+1
    cv = np.asarray(cv)
    count = cv.shape[0]
    if periodic:
        factor, fraction = divmod(count+degree+1, count)
        cv = np.concatenate((cv,) * factor + (cv[:fraction],))
        count = len(cv)
        degree = np.clip(degree,1,degree)
    # If opened, prevent degree from exceeding count-1
    else:
        degree = np.clip(degree,1,count-1)
    # Calculate knot vector
    if periodic:
        kv = np.arange(0-degree,count+degree+degree-1,dtype='int')
    else:
        kv = np.array([0]*degree + list(range(count-degree+1)) + [count-degree]*degree,dtype='int')
    def spl(u, normalize=True):
        """normalize： 終点がu=1になるように正規化する"""
        # Calculate result
        u = np.asarray(u)
        if normalize:
            if periodic:
                u *= count - degree - 1
            else:
                u *= count - degree
        points = np.zeros((len(u), cv.shape[1]))
        for i in range(cv.shape[1]):
            points[:,i] = splev(u, (kv,cv[:,i],degree))
        return points
    return spl


def interpolatedline(p0, p1):
    func = [InterpolatedUnivariateSpline((0,1),(x0,x1),k=1) for (x0,x1) in zip(p0,p1)]
    def line(t):
        return Point([f(t) for f in func])
    return line


def main(xfin=0., yfin=192., zfin=0.5, delta=0., span_fin=50., rchord_fin=80., taper=1., sweep=0.):
    include_fins = True  # 垂直フィンを含める場合True（Falseの場合はクリーン形態）

    wgs_aux = WGS_AUX("Created by STakanashi by makewgs")
    # 主翼のwgsを記述する
    # <editor-fold desc="main_wing">
    np_wing_x = 35
    np_wing_y = 35
    np_wing_z = 4

    chord_wing = 80.
    halfspan_wing = 202.
    naca2410u = partial(naca4digit, 0, 0, 12, surface = 'upper')
    naca2410l = partial(naca4digit, 0, 0, 12, surface = 'lower')
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
    wgs_aux.add_network("wing", wing, 1)

    wingtip = wingtipu.linspace(wingtipl.flip(), np_wing_z)
    wgs_aux.add_network("wingtip", wingtip, 1)
    # </editor-fold>

    # wakeのwgsを記述する
    # <editor-fold desc="wake">
    wake_length = chord_wing * 25.

    wingwake = wing.make_wake(3, wake_length)
    wgs_aux.add_network("wingwake", wingwake, 18)
    # </editor-fold>

    wgs_aux.create_wgs("model3.wgs")

    #aoa = np.arange(-6, 18, 2)
    #aoa = list(zip(*[iter(aoa)] * 3))
    aoa = (2,)
    for (i, a) in enumerate(aoa):
        wgs_aux.create_aux("model3_{}.aux".format(i + 1), a)

    wgs_aux.create_stl("model3.stl")


if __name__ == '__main__':
    main()
