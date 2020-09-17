import os
import shutil
import subprocess
import time
import myutils as ut


def main():
    # for d in ut.iglobm('path/gen_*'):
    #     if os.path.isdir(d):
    #         shutil.rmtree(d)

    # for f in ut.iglobm('result/*.csv'):
    #     os.remove(f)

    # for d in ut.iglobm('result/*'):
    #     if os.path.isdir(d):
    #         shutil.rmtree(d)

    for d in ['path', 'result']:
        if os.path.isdir(d):
            os.rename(d, ut.uniq_path(f'{d}_{ut.strnow("%Y/%m/%d")}'))


    with ut.stopwatch('main') as sw:
        print(time.strftime('%Y/%m/%d %H:%M:%S'))
        result = subprocess.run('dof_new.exe')

        if not result.returncode == 0:
            with ut.EmailIO(None, 'Error: dof_new.exe') as e:
                print(ut.strnow(), file=e)

        elif sw > 3600:
            with ut.EmailIO(None, 'Completed: dof_new.exe') as e:
                print(ut.strnow(), file=e)


if __name__ == '__main__':
    main()
