import glob
from itertools import combinations
from operator import attrgetter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import myutils as ut


def check_dominated(array_s, array_o):
    eq = array_s == array_o
    le = array_s <= array_o
    ge = array_s >= array_o

    if eq.all():
        return 0
    elif le.all():
        return 1
    elif ge.all():
        return -1
    else:
        return 0


def select_front(array):
    array = np.asarray(array)
    n = len(array)

    if n <= 1:
        return array

    front = np.full((n, n), True)
    for i, j in combinations(range(n), 2):
        d = check_dominated(array[i], array[j])
        if d == 1:
            front[j, i] = False
        elif d == -1:
            front[i, j] = False

    return array[front.all(axis=1)]


def hypervolume(array, ref_point):
    if len(array) == 0:
        return 0

    p = None
    s = 0
    array = np.asarray(array)
    array = array[array[:, 0].argsort(), :]
    ref_point = np.asarray(ref_point)

    for a in array:
        if np.any(a > ref_point):
            continue

        elif p is not None:
            w = a[0] - p[0]
            h = ref_point[1] - p[1]
            assert w >= 0
            assert h >= 0
            s += w * h

        p = a

    if p is not None:
        w = ref_point[0] - p[0]
        h = ref_point[1] - p[1]
        assert w >= 0
        assert h >= 0
        s += w * h

    return s


def __test__():
    data = [[1, -1], [3, -2], [4, -4], [10, 10]]
    ref = [5, 0]

    print(hypervolume([], []))

    assert hypervolume(data, ref) == 8


################################################################################

def get_hvs(file):
    def f_():
        for s in steps:
            array = df.query(f'step=={s}')[['obj1', 'obj2']]
            array = np.asarray(array)
            array = select_front(array)
            hv = hypervolume(array, ref_point)
            yield s, hv

    # ref_point = np.array([1.6e7, 9e5])
    # ref_point = np.array([1.5e7, 9e5]) # 1
    # ref_point = np.array([1.5e7, 3e6]) # 2
    ref_point = np.array([1.5e7, 70]) # 3: J, We
    # obj1: 'Cost function, J'
    # obj2: 'Total impulse [Ns]'

    print('read:', file)
    df = pd.read_csv(file)
    steps = sorted(set(df['step']))
    df = df.query("feasible=='T'")
    hvs = np.array(list(f_()))
    return hvs


def plot_hv():
    plt.style.use('ggplot')

    def plot(png):
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci',  axis='x', scilimits=(0, 0))
        ax.ticklabel_format(style='sci',  axis='y', scilimits=(0, 0))
        # plt.show()
        # plt.xlim([0.8e7, 1.6e7])
        # plt.ylim([4e5, 9e5])
        plt.xlabel('Generation')
        plt.ylabel('Hypervolume')
        plt.legend()
        plt.savefig(png, bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

    for mb in ['woMB']:
        hvs = []

        for file in ut.iglobm(f'test4_*_{mb}_*/result/history.csv'):
            label = 'BLX' if 'blx' in file else 'SBX'
            png = f'hypervolume_test4_{mb}_{label}.png'
            hv = get_hvs(file)
            plt.plot(*hv.T, label=label)
            plot(png)
            hvs.append({'label':label, 'hv':hv})

        png = f'hypervolume_test4_{mb}.png'
        for label, hv in hvs.items():
            plt.plot(*hv.T, label=label)
        plot(png)
        # file0 = f'test2_2obj_blx_pop{npop}/result/history.csv'
        # file1 = f'test2_2obj_sbx-fix_pop{npop}/result/history.csv'

        # plt.plot(*get_hvs(file0).T, label='BLX')
        # plt.plot(*get_hvs(file1).T, label='SBX')



if __name__ == '__main__':
    plot_hv()
