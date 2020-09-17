import glob
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd
import myutils as ut


def separate_main(filename):
    def f_(reader):
        columns = None

        for i, df in enumerate(reader):
            print('read:', i, end=' \r')

            if columns is None:
                columns = [s for s in df.columns if not s.startswith('var')]

            step_min = df['step'].min()
            step_max = df['step'].max()

            for n in range(step_min, step_max+1):
                if (n < 1000 and n % 100 == 0) or n % 1000 == 0:
                    df_s = df.query(f'step=={n}')

                    if n in df_table:
                        df_table[n] = pd.concat([df_table[n], df_s],
                                                ignore_index=True)
                    else:
                        df_table[n] = df_s

            for n, df_s in list(df_table.items()):
                if n < step_max:
                    df_s.to_csv(f'{output}/step_{n:05d}.csv')
                    del df_table[n]

            yield df[columns]


    # filename = 'out_so_X_5.csv'
    if os.path.getsize(filename) < 1024**3:
        return
    output = f'{filename[:-4]}_separated'
    if os.path.isdir(output):
        return
    print(os.path.abspath(output))
    os.makedirs(output, exist_ok=True)


    df_table = {}
    reader = pd.read_csv(filename, chunksize=1000)
    df = pd.concat(f_(reader), ignore_index=True)


    for n, df_s in df_table.items():
        df_s.to_csv(f'{output}/step_{n:05d}.csv')
    df.to_csv(f'{output}/result_total.csv')


def plot_result():
    for file in ut.iglobm('**/result_total.csv'):
        png = file[:-4] + '.png'

        if os.path.exists(png):
            os.remove(png)
            # continue

        print('print:', png, end=' \r')

        df = pd.read_csv(file).query("feasible=='T'")
        if len(df) == 0:
            continue

        df.plot(kind='scatter', x='obj2', y='obj1', s=3)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci',  axis='x', scilimits=(0, 0))
        ax.ticklabel_format(style='sci',  axis='y', scilimits=(0, 0))
        # plt.xlim([4e5, 9e5])
        # plt.ylim([0.6e7, 1.5e7])
        plt.xlabel('Total impulse [Ns]')
        plt.ylabel('Cost function, J')
        # plt.legend()
        plt.savefig(png, bbox_inches='tight', pad_inches=0.1)
        plt.close('all')


def main():
    for file in ut.iglobm('**/out_so_X_4.csv'):
        dn = os.path.dirname(file)
        fn = os.path.basename(file)
        with ut.chdir(dn):
            separate_main(fn)


if __name__ == '__main__':
    with ut.chdir(r'G:\2019_1119_着陸飛行経路最適化\landing_opt'):
        plot_result()
