# coding: utf-8
import sys, os, glob, re, subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
import csv


def predict2d(X, Y, Z, lx, ly, lz):
  fig = plt.figure(figsize=(8, 7))
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
  # ax = Axes3D(fig)
  plt.xlabel(lx, size=14)
  plt.ylabel(ly, size=14)
  plt.tick_params(labelsize=14)

  im = plt.pcolor(X, Y, Z, cmap='bwr', linewidth=0)
  cbar = fig.colorbar(im)
  cbar.set_label(lz, size=14)
  # plt.contour(X, Y, Z, zdir='z')
  plt.show()

def predict3d(X, Y, Z, lx, ly, lz):
  fig = plt.figure()
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.7, top=0.95, wspace=0.1, hspace=0.2)
  ax = Axes3D(fig)
  ax.set_xlabel(lx, size=14)
  ax.set_ylabel(ly, size=14)
  ax.set_zlabel(lz, labelpad=30, size=14)
  ax.tick_params(labelsize=14)
  ax.tick_params(axis='z', pad=20)

  # plt.gca().zaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
  # ax.set_aspect(0.2)

  ax.plot_surface(X, Y, Z, cmap='bwr', linewidth=0)
  ax.contour(X, Y, Z, zdir='z', offset=np.min(Z))
  plt.show()
  # ax.imwrite("out.png")

def predict3d2(X, Y, Z, P, Q, R, lx, ly, lz):
  fig = plt.figure()
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.7, top=0.95, wspace=0.1, hspace=0.2)
  ax = Axes3D(fig)
  ax.set_xlabel(lx, size=14)
  ax.set_ylabel(ly, size=14)
  ax.set_zlabel(lz, labelpad=30, size=14)
  ax.tick_params(labelsize=14)
  ax.tick_params(axis='z', pad=20)

  # plt.gca().zaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
  # ax.set_aspect(0.2)

  ax.plot_surface(X, Y, Z, cmap='bwr', linewidth=0)
  ax.plot(P, Q, R, "o")
  ax.contour(X, Y, Z, zdir='z', offset=np.min(Z))
  plt.show()
  # ax.imwrite("out.png")

def predict3d_scatter(X, Y, Z, lx, ly, lz):
  fig = plt.figure()
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.7, top=0.95, wspace=0.1, hspace=0.2)
  ax = Axes3D(fig)
  ax.set_xlabel(lx, size=14)
  ax.set_ylabel(ly, size=14)
  ax.set_zlabel(lz, labelpad=30, size=14)
  ax.tick_params(labelsize=14)
  ax.tick_params(axis='z', pad=20)

  # plt.gca().zaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
  # ax.set_aspect(0.2)

  ax.scatter(X, Y, Z, cmap='bwr', vmin=Z.min(), vmax=Z.max())
  # ax.contour(X, Y, Z, zdir='z', offset=np.min(Z))
  plt.show()
  # ax.imwrite("out.png")

def main(argv):
  if len(argv) == 0:
    print("ファイル番号を指定してください: { 1, 2, 3 }")
    no = sys.stdin.readline().rstrip("\n")
    print("表示形式を指定してください: 2d=>2, 3d=>3, 3d(scatter)=>4")
    dim = sys.stdin.readline()[0]
  elif len(argv) == 1:
    no = argv[1]
    dim = "2"
  else:
    no = argv[0]
    dim = argv[1][0]

  file = "test" + no + ".csv"
  if not os.path.exists(file):
    print("入力ファイルがありません: " + file)
    return

  labels = ["$\\alpha, deg$", "$elv, deg$", "$Ma$", "$x_1$", "$x_2$"]
  x_label = labels[3]
  y_label = labels[4]
  z_label = ["$c_x$", "$c_m$", "$c_z$", "$z$"][int(no) - 1]

  to_grid = dim in "23"

  with open(file, "r") as f:
    X, Y, Z = (np.array(x, dtype=np.float) for x in list(zip(*csv.reader(f)))[0:3])

    if to_grid:
      X, Y, Z = (x.reshape((101, -1)) for x in (X, Y, Z ))

    if dim == "2":
      predict2d(X, Y, Z, x_label, y_label, z_label)
    elif dim == "3":
      predict3d(X, Y, Z, x_label, y_label, z_label)
    elif dim == "4":
      predict3d_scatter(X, Y, Z, x_label, y_label, z_label)

  # with open("table_"+["xl", "ym", "zn"][int(no) - 1]+".csv", "r") as f:
  #   P, Q, R = (np.array(x, dtype=np.float) for x in list(zip(*csv.reader(f)))[0:3])
  # predict3d2(X, Y, Z, P, Q, R, x_label, y_label, z_label)

main(sys.argv[1:])
