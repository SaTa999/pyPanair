import numpy as np
from scipy.interpolate import splev


def naca4digit(m, p, xx, x, c, surface):
    """
    return the coordiantes of a naca 4 digit airfoil
    e.g. for naca2412, m=2, p=4, xx=12
    """
    m /= 100.0
    p /= 10.0
    t = xx / 100.0
    xi = x / c
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # originally -0.1015, but -0.1036 for sharp edge
    yt = 5 * t * (a0 * (xi ** 0.5) + a1 * xi + a2 * (xi ** 2) + a3 * (xi ** 3) + a4 * (xi ** 4))
    if x < 0. or c < x:
        raise ValueError("x = {0} is out of range for c = {1}".format(x, c))
    else:
        pass
    # modify the TE to make the its thickness 0
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
    theta = np.arctan(grad_yc)
    if surface == 'upper':
        pm = 1
    elif surface == 'lower':
        pm = -1
    else:
        raise ValueError('Error! surface must be \'upper\' or \'lower\'')
    xu = c * (xi - pm * yt * np.sin(theta))
    yu = c * (yc + pm * yt * np.cos(theta))
    return [xu, yu]


def bspline(cv, degree=3, periodic=False):
    """ return a function that defines a bezier spline
        cv :      an array of control points
        degree:   degree of the polynominal curve
        periodic: True - closed curve
                  False - open curve
        adopted from http://stackoverflow.com/a/35007804
    """
    # If periodic, extend the point array by count+degree+1
    cv = np.array(cv)
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
        """normalize： normalize the parameter u so, u=1 at the end of the curve"""
        # Calculate result
        u = np.array(u)
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

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

    # cv = np.array([[50., 25.],
    #                [59., 12.],
    #                [50., 10.],
    #                [57., 2.],
    #                [40., 4.],
    #                [40., 14.]])

    n_cpoints = 4
    sep = 1 / (n_cpoints - 2)
    cv = np.zeros((n_cpoints, 2))
    for i in range(1, n_cpoints-1):
        cv[i,0] = np.random.random() * sep + sep * (i-1)
    cv[-1,0] = 1
    cv[1:,1] = np.random.random(cv.shape[0]-1) * 2 -1


    plt.plot(cv[:, 0], cv[:, 1], 'o', label='Control Points')

    us = np.linspace(0, 1, 100)
    for d in range(1, 4):
        res =  bspline(cv, degree=d, periodic=False)(us)
        x = res[:,0]
        y = res[:,1]
        plt.plot(x, y, '-', label='Degree %s' % d)

    plt.minorticks_on()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-1, 1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()