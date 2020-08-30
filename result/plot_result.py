import glob
import os
import time
import traceback
import pandas as pd
from pandas import DataFrame
import numpy as np
from matplotlib import pyplot as plt


def main():
    plt.style.use('ggplot')

    def f_():
        for file in glob.glob('out_*.csv'):
            print('read:', file)
            df = pd.read_csv(file)
            xlim = df['obj1'].min(), df['obj1'].max()
            ylim = df['obj2'].min(), df['obj2'].max()
            yield file, df, np.array([xlim, ylim])


    data = list(f_())
    mn = np.stack([x[2][:, 0] for x in data]).min(axis=0)
    mx = np.stack([x[2][:, 1] for x in data]).max(axis=0)

    for file, df, _ in data:
        print('print:', file)
        png = file[:-4] + '.png'
        # if os.path.exists(png):
        #     continue

        df.plot(kind='scatter', x='obj1', y='obj2')
        plt.xlim([0, mx[0]*1.1])
        plt.ylim([0, mx[1]*1.1])
        plt.xlabel('Total impulse [Ns]')
        plt.ylabel('Cost function, J')
        plt.savefig(png)
        plt.close('all')


if __name__ == '__main__':
    try:
        main()
    except:
        traceback.print_exc()
    finally:
        time.sleep(2)
