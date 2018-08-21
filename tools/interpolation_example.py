import numpy as np

import h5py

# helper to open a dataset
from simdataset import simdataset

# to handle path
import pathlib

# to handle arguments
import argparse
import re

from math import cos, sin, sqrt, pi, pow
from utilities import spherical
from icogrid import ico
import time

# argument parsing
parser = argparse.ArgumentParser(
    description='THOR dataset interpolation example')

parser.add_argument("-f",
                    "--folder",
                    action='store',
                    metavar='FOLDER',
                    help='folder to read data from',
                    default="../results/")

args = parser.parse_args()

# get dataset folder
folder = pathlib.Path(args.folder)

# get data
dataset = simdataset(folder)
planets = dataset.get_planets()
if len(planets) < 1:
    print("no planet found in dataset")
    exit(-1)

dataset.select_planet(planets[0])
dataset.select_scalar_data_type("Momentum_norm")
dataset.select_vector_data_type("Momentum")


num_samples = dataset.get_num_datasets()
print("num datasets:", num_samples)
num_points, lvl = dataset.get_dim()


# get grid points (array of lon lat coordinates)
# this returnes a reshaped array from the h5py array, for easier indexing:
# np.array(lonlat).reshape(num_points, 2)
lonlat = dataset.get_grid()

# get some data for first file, reshaped as (num_points, num_levels, num:Data_dim)
data = dataset.dataloader.get_data(0, "Rho")

# compute interpolation points
# we'll interpolate on points along a great circle defined by its
# perpendicular
num_interp = 360


interp_vectors = np.zeros((num_interp, 3))

# choose a vector perpendicular to a great circle, dephined, by theta, phi
theta = 25.0*pi/(2.0*90.0)
phi = 10.0*pi/(2.0*90.0)
r = spherical(1.0, theta, phi)
# compute perpendicular in z=0 plane
p_z = 0.0
p_y = sqrt(1.0/(1+pow(r[1]/r[0], 2.0)))
p_x = -r[1]/r[0]*p_y

# compute two vectors in plane of great circle
p1 = np.array((p_x, p_y, p_z))
p2 = np.cross(p1, r)

# compute all the interpolation points on the great circle
for i in range(num_interp):
    alpha = (i/(num_interp))*pi*2.0
    interp_vectors[i, :] = cos(alpha)*p1+sin(alpha)*p2


g = int(pow((num_points - 2)/10, 1/4)) - 2

print("Start computation for barycentric coordinates")
start = time.time()

# prepare the grid to get barycentric coordinates
icogrid = ico(g, lonlat)

# compute the coordinates on the mesh for each ray
dot, t_i, u_i, v_i = icogrid.barycentric_coordinates.get_barycentric_coordinates(
    interp_vectors)

# compute u/v coordinates
u = u_i/dot
v = v_i/dot

# get the triangles that intersected with rays in the direction of the
# ray (and not behind), and with intersection inside the triangle,
# i.e. 0 < u < 1, 0 < v < 1, u+v < 1
condition = np.logical_and.reduce((dot > 0.0,
                                   u > 0.0,
                                   u < 1.0,
                                   v > 0.0,
                                   v < 1.0,
                                   u+v < 1.0))

# geta numpy masc for the interesting points

w = np.where(condition)
stop = time.time()
print("Barycentric coordinates done: {} s".format(stop-start))

level = 0

# loop on all intersections:
for i in range(w[0].shape[0]):
    # get triangle index in grid
    triangle_idx = w[0][i]
    # get index of triangle edges in lonlat grid
    triangle = icogrid.triangles[triangle_idx]
    i1 = triangle[0]
    i2 = triangle[1]
    i3 = triangle[2]

    # compute vertices of triangle from lonlat grid
    p1 = spherical(1.0, lonlat[i1, 1], lonlat[i1, 0])
    p2 = spherical(1.0, lonlat[i2, 1], lonlat[i2, 0])
    p3 = spherical(1.0, lonlat[i3, 1], lonlat[i3, 0])

    # get u, v values for this intersection
    u_p = u[w[0][i], w[1][i]]
    v_p = v[w[0][i], w[1][i]]

    # compute triangle edges
    v0 = p2 - p1
    v1 = p3 - p1

    # compute intersection point
    pk = p1 + u_p*v0 + v_p * v1

    # get matching data from array for each corner
    data1 = data[i1, level]
    data2 = data[i2, level]
    data3 = data[i3, level]

    # data can be interpolated from u,v,w weights
    w_p = 1-u-v

    # P = w*p1 + u*p2 + v*p3
    d = w_p*data1 + u_p * data2 + v_p * data3
