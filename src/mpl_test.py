import ctypes
import random
from itertools import chain, islice
from operator import attrgetter

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


################################################################################

def get_f_init(cdll):
    ''' ラッパー1: init
    '''
    func_ptr = getattr(cdll, 'init')
    func_ptr.argtypes = []
    func_ptr.restype = ctypes.c_void_p
    def f_():
        func_ptr()
    return f_


def get_f_calc(cdll, n_data=2):
    ''' ラッパー2: calc
    '''
    func_ptr = getattr(cdll, 'calc')
    func_ptr.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        ctypes.c_int32,
        ctypes.c_double,
        np.ctypeslib.ndpointer(dtype=np.float64)
    ]
    func_ptr.restype = ctypes.c_void_p
    def f_(args, d):
        if not isinstance(args, np.ndarray):
            args = np.array(args)
        n = ctypes.c_int32(len(args))
        d = ctypes.c_double(d)
        result = np.zeros((len(args), n_data), dtype=np.float64)
        func_ptr(args, n, d, result)
        return result
    return f_


################################################################################

class Point2d(object):

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
        self.pinned = False

    def __getitem__(self, key):
        return (self.x, self.y)[key]

    def __lt__(self, other):
        return self.x < other.x

    @property
    def array(self):
        return np.array([self.x, self.y])

    def setloc(self, x, y):
        self.x = 0 if self.pinned else max(x, 0)
        self.y = y
        return self


################################################################################

def rand_it():
    while True:
        yield random.random()


def onmousemove(event):
    print('moved', event.x, event.y, event.xdata, event.ydata, end='\r')


def get_points(n, x_max=1):
    x = np.linspace(0, x_max, n)
    # y = np.sin(2 * np.pi * x)
    y = np.zeros_like(x)
    points = list(map(Point2d, x, y))
    points[0].pinned = True
    return points

    # return list(islice(map(Point2d, rand_it(), rand_it()), n))


################################################################################

class Handler(object):
    def __init__(self, data, func, fig, axes):
        self.data = data
        self.func = func
        self.fig = fig
        self.axes = axes
        self.picked = None

        fig.canvas.mpl_connect('pick_event', self.onpick)
        fig.canvas.mpl_connect('motion_notify_event', self.onmove)
        fig.canvas.mpl_connect('button_release_event', self.onrelease)
        fig.canvas.mpl_connect('button_press_event', self.onclick)

        self.calc_data()
        self.calc_line()
        self.plot_data()
        self.plot_result()
        # self.fig.legend()
        plt.show()

    def onpick(self, event):
        if self.picked is not None:
            return

        # artist = event.artist
        # xdata = artist.get_xdata()
        # ydata = artist.get_ydata()
        ind = event.ind[0]
        self.picked = self.data[ind]
        print('picked:', ind, end=' '*20+'\r')

    def onmove(self, event):
        if self.picked is None:
            return

        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            return
        print('moving', end=' '*20+'\r')

        self.picked.setloc(x, y)
        self.calc_data()
        self.plot_data()

    def onrelease(self, event):
        print('released', end=' '*20+'\r')
        self.calc_line()
        self.plot_data()
        self.picked = None

    def onclick(self, event):
        ''' 右クリック '''

        if not event.button == 3:
            return
        print('clicked', end=' '*20+'\r')
        self.plot_result()

    def calc_data(self):
        self.data.sort()
        array = np.array(list(map(attrgetter('array'), self.data))).T
        self.array = array

    def plot_data(self):
        self.axes[-1].cla()

        # x_margin = 0.1 * (self.x_max - self.x_min)
        # y_margin = 0.1 * (self.y_max - self.y_min)
        # self.ax.set_xlim((self.x_min - x_margin, self.x_max + x_margin))
        # self.ax.set_ylim((self.y_min - y_margin, self.y_max + y_margin))
        self.axes[-1].set_ylim((-40, 40))
        self.axes[-1].set_xlabel('time [t]')
        self.axes[-1].set_ylabel('elevator [°]')

        self.axes[-1].plot(*self.array, 'o', label='Raw data', picker=5)
        self.axes[-1].plot(self.linex, self.liney, '-', label='Cubic spline')
        self.fig.canvas.draw()

    def calc_line(self, n=100):
        spline = interp1d(*self.array, kind='cubic')
        linex = np.linspace(self.x_min, self.x_max, n)
        liney = spline(linex)
        self.linex = linex
        self.liney = liney

    def get_linedata(self, dx, x_max=float('inf')):
        spline = interp1d(*self.array, kind='cubic')
        linex = np.arange(0, min(self.x_max, x_max), dx)
        liney = spline(linex)
        return liney

    def plot_result(self):
        x, z, a, t = self.func(self.get_linedata)

        self.axes[0].cla()
        self.axes[0].plot(x, z, '-')
        self.axes[0].set_xlabel('x [m]')
        self.axes[0].set_ylabel('z [m]')

        self.axes[1].cla()
        self.axes[1].plot(a, t, '-')
        self.axes[1].set_xlabel('time [s]')
        self.axes[1].set_ylabel('angle of attack [°]')
        self.fig.canvas.draw()

    @property
    def x_min(self):
        return self.array[0, :].min()

    @property
    def x_max(self):
        return self.array[0, :].max()

    @property
    def y_min(self):
        return self.array[1, :].min()

    @property
    def y_max(self):
        return self.array[1, :].max()


def plot_points_main():

    # CDLLインスタンス作成
    libname = 'dof.dll'
    loader_path = '.'
    cdll = np.ctypeslib.load_library(libname, loader_path)

    # 関数取得
    f_init = get_f_init(cdll)
    f_calc = get_f_calc(cdll, n_data=6)

    # 初期化: 開始時刻設定
    f_init()

    def func(callback):
        dt = 0.01
        elv = callback(dt)
        orbit = f_calc(elv, dt)
        times = np.arange(len(elv)) * dt
        # return orbit[:, 0], orbit[:, 1], times, orbit[:, 5]
        return times, orbit[:, 4], times, orbit[:, 5]


    fig, axes = plt.subplots(3, 1, figsize=(10, 7))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=1,
                        wspace=0, hspace=0.25)

    points = get_points(10, x_max=5)
    handler = Handler(points, func=func, fig=fig, axes=axes)

    # plot_points_with_spline(points, ax)
    # fig.legend()


    # fig.canvas.mpl_connect('button_press_event', onclick)
    # fig.canvas.mpl_connect('pick_event', onpick)
    # fig.canvas.mpl_connect('motion_notify_event', onmove)
    # fig.canvas.mpl_connect('button_release_event', onrelease)
    # plt.show()


################################################################################

def __test__():
    def rand_it():
        while True:
            yield random.random()

    point_it = map(Point2d, rand_it(), rand_it())
    points = sorted(islice(point_it, 10))
    return np.array([[0, 0]] + [[p.x, p.y] for p in points] + [[1, 1]])


def __test__():
    fig, ax = plt.subplots()
    # t = np.linspace(0, 1, 41)
    # x = 2 * np.pi * t
    # y = np.sin(x ** 2)
    data = __test__()
    x, y = data[:, 0], data[:, 1]

    spline = interp1d(x, y, kind='cubic')

    t_ = np.linspace(0, 1, 1000)
    # x_ = 2 * np.pi * t_
    # y_ = spline(x_)
    # yt = np.sin(x_ ** 2)
    x_ = t_
    y_ = spline(x_)

    ax.plot(x, y, 'o', label='Raw data')
    ax.plot(x_, y_, '-', label='Cubic spline')
    # ax.plot(x_, yt, '-', label='True line')
    fig.legend()

    fig.canvas.mpl_connect('motion_notify_event', onmousemove)
    plt.show()


def __test__():
    # CDLLインスタンス作成
    libname = 'dof.dll'
    loader_path = '.'
    cdll = np.ctypeslib.load_library(libname, loader_path)

    # 関数取得
    f_init = get_f_init(cdll)
    f_calc = get_f_calc(cdll)

    # 初期化: 開始時刻設定
    f_init()

    elv = np.ones(100, dtype=np.float64)
    dt = 0.01

    results = f_calc(elv, dt)
    print(results)

    plt.plot(*results.T)
    plt.show()


def main():
    # return __test__()
    plot_points_main()


if __name__ == '__main__':
    main()
