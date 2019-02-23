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


class Handler(object):
    def __init__(self, points, func, fig, axes):
        self.points = points
        self.func = func
        self.fig = fig
        self.axes = axes
        self.picked = None
        self.busy = False

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

    def __len__(self):
        return len(self.points)

    def onpick(self, event):
        ''' 制御点クリック
        '''
        if self.picked is not None:
            return

        # artist = event.artist
        # xdata = artist.get_xdata()
        # ydata = artist.get_ydata()

        ind = event.ind[0]
        self.picked = self.points[ind]
        print('picked:', ind, end=' '*20+'\r')

        if event.mouseevent.button == 1:
            # 左クリック
            pass

        if event.mouseevent.button == 2:
            # ホイールクリック
            if ind == len(self) - 1:
                self.picked.x *= 1.2
                self.calc_data()
                self.plot_data()

        elif event.mouseevent.button == 3:
            # 右クリック
            if ind and len(self) > 1:
                self.points.pop(ind)
                self.calc_data()
                self.plot_data()

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
        if self.picked is None:
            return

        print('released', end=' '*20+'\r')
        self.picked = None
        self.calc_line()
        self.plot_data()
        self.plot_result() # 経路計算

    def onclick(self, event):
        if self.picked is not None:
            return

        if event.button == 3:
            # 右クリック
            x = event.xdata
            y = event.ydata
            if x is None or y is None:
                return

            print('clicked', end=' '*20+'\r')
            point = Point2d(x, y)
            self.points.append(point)
            self.picked = point
            self.calc_data()
            self.plot_data()

    def calc_data(self):
        self.points.sort()
        self.array = np.array(list(map(attrgetter('array'), self.points))).T

    def plot_data(self):
        self.axes[-1].cla()
        # x_margin = 0.1 * (self.x_max - self.x_min)
        # y_margin = 0.1 * (self.y_max - self.y_min)
        # self.ax.set_xlim((self.x_min - x_margin, self.x_max + x_margin))
        # self.ax.set_ylim((self.y_min - y_margin, self.y_max + y_margin))
        self.axes[-1].set_ylim((-40, 40))
        self.axes[-1].set_xlabel('time [s]')
        self.axes[-1].set_ylabel('elevator [°]')
        self.axes[-1].plot(*self.array, 'o', label='Raw data', picker=5)
        self.axes[-1].plot(self.linex, self.liney, '-', label='Cubic spline')
        self.axes[-1].patch.set_alpha(0)
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
        if self.busy:
            return

        try:
            self.busy = True
            for i, data in enumerate(self.func(self.get_linedata)):
                self.axes[i].cla()
                self.axes[i].plot(data['x'], data['y'], '-')
                self.axes[i].set_xlabel(data.get('lx'))
                self.axes[i].set_ylabel(data.get('ly'))
                self.axes[i].patch.set_alpha(0)
            self.fig.canvas.draw()

        finally:
            self.busy = False

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


################################################################################

def get_points(n, x_max=1):
    ''' 入力の初期値を生成
    '''
    x = np.linspace(0, x_max, n)
    y = np.zeros_like(x) - 10
    points = list(map(Point2d, x, y))
    points[0].pinned = True
    return points


def plot_points_main():
    ''' CDLLインスタンス作成
    '''
    libname = 'libdof.dll'
    loader_path = '.'
    cdll = np.ctypeslib.load_library(libname, loader_path)

    # 関数取得
    f_init = get_f_init(cdll)
    f_calc = get_f_calc(cdll, n_data=7)

    # 初期化: 開始時刻設定
    f_init()

    def func(callback):
        dt = 0.01
        elv = callback(dt)
        orbit = f_calc(elv, dt)
        times = np.arange(len(elv)) * dt
        # x0, y0 = orbit[:, 0], orbit[:, 1]
        # data0 = {'x': times, 'y': orbit[:, 5], 'lx': 'x [m]', 'ly': 'z [m]'}
        data = (
            {'x': orbit[:, 0], 'y': orbit[:, 1], 'lx': 'x [m]', 'ly': 'z [m]'},
            {'x': times, 'y': orbit[:, 5], 'lx': 'time [t]', 'ly': '$\\alpha$ [°]'},
            {'x': times, 'y': orbit[:, 4], 'lx': 'time [t]', 'ly': '$\\theta$ [°]'},
            # {'x': times, 'y': orbit[:, 2], 'lx': 'time [t]', 'ly': '$c_m$ [-]'},
            {'x': times, 'y': orbit[:, 6], 'lx': 'time [t]', 'ly': 'Ma [-]'}
        )

            # self.axes[1].cla()
            # self.axes[1].plot(a, t, '-')
            # self.axes[1].set_xlabel('time [s]')
            # # self.axes[1].set_ylabel('angle of attack [°]')
            # # self.axes[1].set_ylabel('$\\theta$ [°]')
            # self.axes[1].set_ylabel('value')
        return data
        # return times, orbit[:, 4], times, orbit[:, 5] # aoa, cx


    fig, axes = plt.subplots(5, 1, figsize=(7, 8))
    fig.subplots_adjust(left=0.12, bottom=0.06, right=0.95, top=1,
                        wspace=0, hspace=0.4)
    fig.patch.set_alpha(0.5)

    points = get_points(10, x_max=100)
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
