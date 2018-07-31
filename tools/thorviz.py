import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import os
import h5py

import pathlib
import argparse
import re

from math import sqrt

from PyQt5 import QtGui, QtCore, uic
from PyQt5.QtWidgets import QMainWindow, QApplication, qApp

from PyQt5.QtCore import QByteArray, QIODevice, Qt, QTimer, pyqtSignal, pyqtSlot
from simdataset import simdataset


# argument parsing
parser = argparse.ArgumentParser(description='THOR dataset visualisation tool')

parser.add_argument("-f",
                    "--folder",
                    action='store',
                    metavar='FOLDER',
                    help='folder to read data from',
                    default="../results/")

args = parser.parse_args()

folder = pathlib.Path(args.folder)

#grid_dataset = "../results/grid_ref/grid_test_6_grid.h5"

#print("grid def")
#grid = h5py.File(grid_dataset)
# for k, v in grid.items():
#    print(k, v)

#xyz = grid['xyz']
#maps = grid['maps']

#num_points = int(len(xyz)/3)

dataset = simdataset(folder)
planets = dataset.get_planets()
if len(planets) < 1:
    print("no planet found in dataset")
    exit(-1)

dataset.select_planet(planets[0])
dataset.select_scalar_data_type("Mh")
dataset.select_vector_data_type("Mh")


class ThorVizWindow(QMainWindow):
    def __init__(self):
        super(ThorVizWindow, self).__init__()
        uic.loadUi('VizWin.ui', self)

        num_samples = dataset.get_num_datasets()
        print("num datasets:", num_samples)
        pts, lvl = dataset.get_dim()

        self.animation_slider.setMinimum(0)
        self.animation_slider.setMaximum(num_samples - 1)
        self.animation_slider.setSingleStep(1)
        self.animation_spinbox.setMaximum(num_samples - 1)

        self.level_slider.setMinimum(0)
        self.level_slider.setMaximum(lvl - 1)
        self.level_slider.setSingleStep(1)
        self.level_spinbox.setMaximum(lvl - 1)
        self.vizGL.set_grid_data(dataset)
        #self.vizGL.set_grid(grd, neighbours, colors, moments)

        self.actionQuit.triggered.connect(qApp.quit)

        self.actionZero_Position.triggered.connect(self.vizGL.reset_position)
        self.show()

        # self.animation_slider.setMaximum(colors.shape[0])

    def set_image(self, i):
        self.vizGL.set_image(i)

    def show_data(self, b):
        self.animation_slider.setEnabled(b)


app = QApplication([])
window = ThorVizWindow()


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):

        sys.exit(app.exec_())
