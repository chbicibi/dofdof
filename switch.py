import os
import time
import myutils as ut


def get_progress(path):
    if not os.path.isdir(path):
        raise FileNotFoundError(path)

    with ut.chdir(path):
        if not os.path.exists('progress.txt'):
            return 0

        with open('progress.txt', 'r') as f:
            return int(f.read().strip())


def switch(path, i):
    if not os.path.isdir(path):
        raise FileNotFoundError(path)

    with ut.chdir(path):
        if i == 0:
            if os.path.exists('pause'):
                os.remove('pause')

        else:
            if not os.path.exists('pause'):
                with open('pause', 'w'):
                    pass


def main():
    path_table = {
        'test3_2obj_blx_pop30': 0,
        'test3_2obj_sbx_pop30': 0
    }

    try:
        while True:
            for path in list(path_table):
                path_table[path] = get_progress(path)
                if path_table[path] < 0:
                    del path_table[path]

            if not path_table:
                break

            path_list = [p for p, n in sorted(path_table.items(),
                                              key=lambda a: a[1])]
            for i, path in enumerate(path_list):
                switch(path, i)
                if i == 0:
                    print(path, path_table[path], ' '*20, end='\r')

            time.sleep(300)

    finally:
        for path in path_table:
            switch(path, 0)


if __name__ == '__main__':
    main()
