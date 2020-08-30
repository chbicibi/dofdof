import argparse
import glob
import os
import time
import traceback
import pandas as pd
from pandas import DataFrame
import numpy as np
from matplotlib import pyplot as plt


def test():
    file = glob.glob('flight_*.csv')[0]
    df = pd.read_csv(file)
    df.plot(x='time', y='thrust')
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', action='store_true',
                        help='test')
    args = parser.parse_args()
    # plt.style.use('ggplot')

    if args.test:
        test()
        return

    def f_():
        for file in glob.glob('flight_*.csv'):
            print('read:', file)
            df = pd.read_csv(file)
            xlim = df['x'].min(), df['x'].max()
            ylim = df['z'].min(), df['z'].max()
            yield file, df, np.array([xlim, ylim])

    data = list(f_())
    mn = np.stack([x[2][:, 0] for x in data]).min(axis=0)
    mx = np.stack([x[2][:, 1] for x in data]).max(axis=0)

    # for file in glob.glob('flight_*.csv'):
    for file, df, _ in data:
        print('print:', file)
        png = file[:-4] + '.png'
        # if os.path.exists(png):
        #     continue

        # df = pd.read_csv(file)

        df.plot(x='x', y='z')
        plt.xlim([0, mx[0]*1.1])
        plt.ylim([0, mx[1]*1.1])
        plt.xlabel('Distance [m]')
        plt.ylabel('Altitude [m]')
        plt.savefig(png)
        plt.close('all')


if __name__ == '__main__':
    try:
        main()
    except:
        traceback.print_exc()
    finally:
        time.sleep(2)
