import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import os
import h5py

import pathlib
import argparse
import re

from PyQt5 import QtGui, QtCore, uic
from PyQt5.QtWidgets import QMainWindow, QApplication

from PyQt5.QtCore import QByteArray, QIODevice, Qt, QTimer, pyqtSignal, pyqtSlot

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


# read and sort all files
datasets = folder.glob('esp_output_*_*.h5')

dataset_re = re.compile('esp_output_(.*)_(\d+)')


planets = {}

for f in datasets:
    match = dataset_re.match(f.stem)
    if match is not None:
        basename = match.group(1)
        number = match.group(2)
        if basename not in planets:
            planets[basename] = {'datasets': [(number, f)]}
        else:
            planets[basename]['datasets'].append((number, f))

for planet, values in planets.items():
    # sort them numericaly
    d = values['datasets']
    values['datasets'] = [a[1] for a in sorted(d, key=lambda k:int(k[0]))]

    # open grid file
    grid = folder / ('esp_output_grid_{}.h5'.format(planet))
    planets[planet]['grid'] = grid
    # open plaet file
    planet_def = folder / ('esp_output_{}.h5'.format(planet))
    planets[planet]['def'] = planet_def


grd = None
neighbours = None
for planet, files in planets.items():
    print("Planet: ", planet)
    print("Number of datafiles:", len(files['datasets']))
    num_samples = len(files['datasets'])
    print("Planet def")
    planet_def = h5py.File(files['def'])
    for k, v in planet_def.items():
        print(k, v[0])

    print("grid def")
    grid = h5py.File(files['grid'])
    for k, v in grid.items():
        print(k, v)

    num_levels = int(grid['nv'][0]) + 8
    num_datas = 0
    if len(files['datasets']) > 0:
        print("datafiles def")
        datafile0 = h5py.File(files['datasets'][0])
        num_datas = len(datafile0['Pressure'])

        for k, v in datafile0.items():
            print(k, v)

    if grd is None:
        grd = grid['lonlat']
        neighbours = grid['pntloc']
    num_points = len(grd)//2
    ground_moment = np.zeros(
        (len(files['datasets']), num_datas//(num_levels), 3), dtype=np.float32)
    print(ground_moment.shape)
    for i in range(len(files['datasets'])):
        d = h5py.File(files['datasets'][i])
        ground_moment[i, :, :] = np.array(d['Mh']).reshape(
            num_points, num_levels, 3)[:, 0, :]

    pressure = np.zeros(
        (len(files['datasets']), num_datas//(num_levels)), dtype=np.float32)

    for i in range(len(files['datasets'])):
        d = h5py.File(files['datasets'][i])
        pressure[i, :] = np.array(d['Pressure']).reshape(
            num_points, num_levels)[:, 0]
    print(pressure.shape)
#    for n in neighbours:
#        print(n)
#    print(neighbours)

colors = np.zeros((num_samples, len(grid['lonlat'])//2, 3), dtype=np.float32)
# colors[:, :, 0] = pressure/np.max(pressure)
k = 0.2
colors[:, :, 0] = k + (1.0-k)*ground_moment[:, :, 0] / \
    np.max(ground_moment[:, :, 0])
colors[:, :, 1] = k + (1.0-k)*ground_moment[:, :, 1] / \
    np.max(ground_moment[:, :, 1])
colors[:, :, 2] = k + (1.0-k)*ground_moment[:, :, 2] / \
    np.max(ground_moment[:, :, 2])
# for i in range(n):
#     colors[i, :, 0] = i/n

print(colors)


moments = np.zeros((num_samples, len(grid['lonlat'])//2, 3), dtype=np.float32)
# colors[:, :, 0] = pressure/np.max(pressure)
k = 0.2
moments[:, :, :] = ground_moment[:, :, :]


class ThorVizWindow(QMainWindow):
    def __init__(self):
        super(ThorVizWindow, self).__init__()
        uic.loadUi('VizWin.ui', self)

        self.vizGL.set_grid(grd, neighbours, colors, moments)
        self.show()

        self.animation_slider.setMaximum(colors.shape[0])

    def set_image(self, i):
        pass


app = QApplication([])
window = ThorVizWindow()


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):

        sys.exit(app.exec_())
