import glob
import os
import re
import time
import traceback
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import myutils as ut


def plot_all():
    plt.style.use('ggplot')

    def f_():
        for d in ['test4_*_woMB_*']:
            for file in glob.glob(f'{d}/result/history.csv'):
                print('read:', file, end=' \r')
                df = pd.read_csv(file)#.query("step<=1000")
                # xlim = df['obj1'].min(), df['obj1'].max()
                # ylim = df['obj2'].min(), df['obj2'].max()
                yield file, df


    # data = list(f_())
    # mn = np.stack([x[2][:, 0] for x in data]).min(axis=0)
    # mx = np.stack([x[2][:, 1] for x in data]).max(axis=0)
    # data = f_()

    for file, df in f_():
        # png = file.split('/')[0] + '_all.png'
        png = re.search(r'.+?(?=\\|/|$])', file)[0] + '_all.png'
        print('print:', png, end=' \r')

        with ut.chdir('.'):
            if os.path.exists(png):
                os.remove(png)

            ax = plt.gca()

            # ds = df.query("feasible=='F'")
            # if len(ds) > 0:
            #     ds.plot(kind='scatter', x='obj1', y='obj2', c='red',
            #             label='Infeasible', ax=ax)

            label = 'BLX' if 'blx' in file else 'SBX'

            ds = df.query("feasible=='T'")
            if len(ds) > 0:
                ds.plot(kind='scatter', x='obj2', y='obj1', s=2,# c='blue',
                        label=label, ax=ax)

            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # ax.ticklabel_format(style='sci',  axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # plt.xlim([4e5, 9e5])
            # plt.ylim([0.6e7, 1.5e7])
            # plt.xlabel('Total impulse [Ns]')
            plt.xlabel('$W_e$ [m/s]')
            plt.ylabel('Cost function, J')
            plt.legend()
            plt.savefig(png, bbox_inches='tight', pad_inches=0.1, dpi=150)
            plt.close('all')


def plot_obj():
    def f_():
        for file in glob.glob('out_*.csv'):
            print('read:', file, end=' \r')
            df = pd.read_csv(file)
            xlim = df['obj1'].min(), df['obj1'].max()
            ylim = df['obj2'].min(), df['obj2'].max()
            yield file, df, np.array([xlim, ylim])

    data = list(f_())
    mn = np.stack([x[2][:, 0] for x in data]).min(axis=0)
    mx = np.stack([x[2][:, 1] for x in data]).max(axis=0)
    # data = f_()

    for file, df, _ in data:
        png = file[:-4] + '.png'
        print('print:', png, end=' \r')

        with ut.chdir('obj_'):
            if os.path.exists(png):
                os.remove(png)

            ax = plt.gca()

            ds = df.query("feasible=='F'")
            if len(ds) > 0:
                ds.plot(kind='scatter', x='obj2', y='obj1', c='red', s=3,
                        label='Infeasible', ax=ax)

            ds = df.query("feasible=='T'")
            if len(ds) > 0:
                ds.plot(kind='scatter', x='obj2', y='obj1', c='blue', s=3,
                        label='Feasible', ax=ax)

            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.xlim([0, mx[1]*1.1])
            plt.ylim([0, mx[0]*1.1])
            # plt.xlabel('Total impulse [Ns]')
            plt.xlabel('$W_e$ [m/s]')
            plt.ylabel('Cost function, J')
            plt.legend()
            plt.savefig(png, bbox_inches='tight', pad_inches=0.1)
            plt.close('all')


def plot_each():
    plt.style.use('ggplot')

    # for d in ['test3_2obj_blx_pop30', 'test3_2obj_sbx_pop30']:
    for d in ut.iglobm('test4_*'):
        print(d)
        with ut.chdir(f'{d}/result'):
            plot_obj()


def plot_grad():
    plt.style.use('ggplot')

    file = 'result.csv'
    df = pd.read_csv(file)

    if 'wMB' in os.path.abspath(file):
        label = 'w/ MB'
    else:
        label = 'w/o MB'

    df.query("step<=500").plot(x='step', y='J', label=label)
    plt.ylim([0, 1.7e6])
    plt.xlabel('step')
    plt.ylabel('Cost function, J')
    plt.legend()
    plt.savefig('result.png', bbox_inches='tight', pad_inches=0.1, dpi=150)


def main():
    # plot_all()
    # plot_each()
    plot_grad()


if __name__ == '__main__':
    try:
        main()
    except:
        traceback.print_exc()
    finally:
        time.sleep(2)
