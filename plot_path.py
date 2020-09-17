import argparse
import glob
import os
import time
import traceback
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import myutils as ut


def test():
    file = glob.glob('gen_00600/flight_*.csv')[0]
    df = pd.read_csv(file)
    # df.plot(x='time', y='thrust')

    fig, axes = plt.subplots(3)
    df.plot(x='x', y='z', ax=axes[0])
    df.plot(x='x', y='elv', ax=axes[1])
    df.plot(x='x', y='thrust', ax=axes[2])

    plt.show()


def plot():
    def f_():
        for file in glob.glob('flight_*.csv'):
            print('read:', file, end=' \r')
            df = pd.read_csv(file)
            lims = [
                [df['x'].min(), df['x'].max()],
                [df['z'].min(), df['z'].max()],
                [df['time'].min(), df['time'].max()],
                [df['elv'].min(), df['elv'].max()],
                [df['we'].min(), df['we'].max()]]
            yield file, df, np.array(lims)

    data = list(f_())
    mn = np.stack([x[2][:, 0] for x in data]).min(axis=0)
    mx = np.stack([x[2][:, 1] for x in data]).max(axis=0)
    up = (mx + mn) / 2 + np.maximum((mx - mn) * 0.55, 1)
    lw = (mx + mn) / 2 - np.maximum((mx - mn) * 0.55, 1)

    for file, df, _ in data:
        png = file[:-4] + '.png'
        print('print:', png, end=' \r')
        if os.path.exists(png):
            os.remove(png)

        fig, axes = plt.subplots(3)

        ax = axes[0]
        df.plot(x='x', y='z', ax=ax)
        ax.set_xlim([0, mx[0]*1.1])
        ax.set_ylim([0, mx[1]*1.1])
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Altitude [m]')

        ax = axes[1]
        df.plot(x='x', y='we', ax=ax)
        ax.set_xlim([0, mx[0]*1.1])
        ax.set_ylim([lw[4], up[4]])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('$W_e$ [m/s]')

        ax = axes[2]
        df.plot(x='x', y='elv', ax=ax)
        ax.set_xlim([0, mx[0]*1.1])
        ax.set_ylim([lw[3], up[3]])
        ax.set_xlabel('Distance [m]')
        # ax.set_xlabel('Time [s]')
        ax.set_ylabel('Elevator [deg]')

        # df.plot(x='x', y='thrust', ax=axes[2])
        # axes[2].set_xlim([0, mx[0]*1.1])
        # axes[2].set_ylim([lw[4], up[4]])
        # axes[2].set_xlabel('Distance [m]')
        # # axes[2].set_xlabel('Time [s]')
        # axes[2].set_ylabel('Thrust [N]')

        plt.savefig(png, bbox_inches='tight', pad_inches=0.1)
        plt.close('all')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', action='store_true',
                        help='test')
    args = parser.parse_args()
    # plt.style.use('ggplot')

    if args.test:
        test()
        return

    for d in glob.glob('test4_2obj_*_pop30/path/gen_*'):
        if os.path.isdir(d):
            print(d)
            with ut.chdir(d):
                # if glob.glob('*.png'):
                #     continue
                plot()


if __name__ == '__main__':
    try:
        main()
    except:
        traceback.print_exc()
    finally:
        time.sleep(2)
