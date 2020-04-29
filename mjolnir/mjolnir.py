#!/usr/bin/env python3
import numpy as np
import argparse
from importlib import reload
import time
import mjolnir_plot_helper as mph

first = time.time()
###########################################################################
#
# Options
#
# pview   plot
# uver    Averaged (time and longitude) zonal winds.
#         2D Map Latitude Vs Pressure.
#
# Tver    Averaged (time and longitude) temperatures.
#         2D Map Latitude Vs Pressure.
#
# Tulev   Averaged (time) temperature and wind field.
#         2D Map Longitude Vs Latitude.
#
# PTver   Averaged (time and longitude) potential temperatures.
#         2D Map Latitude Vs Pressure.
#
# ulev    Averaged (time) wind fields in long and lat.
#         2D Map Longitude Vs Latitude.
#
# PVver   Averaged (time and longitude) potential vorticity.
#         2D Map Latitude Vs Pressure.
#
# PVlev   Averaged (time) potential vorticity.
#         2D Map Longitude Vs Latitude.
###########################################################################

parser = argparse.ArgumentParser()
parser.add_argument('pview', metavar='nview', nargs='*', help='Type of plot to make')
parser.add_argument("-f", "--file", nargs=1, default=['results'], help='Results folder to use for plotting')
parser.add_argument("-s", "--simulation_ID", nargs=1, default=['auto'], help='Name of simulation (e.g., planet name)')
parser.add_argument("-i", "--initial_file", nargs=1, default=[10], type=int, help='Initial file id number (integer)')
parser.add_argument("-l", "--last_file", nargs=1, default=['init'], type=int, help='Last file id number (integer)')
parser.add_argument("-lev", "--horizontal_lev", nargs=1, default=[2.5e2], help='Horizonal level to plot in temperature/velocity/vorticity field (mbar or km)')
parser.add_argument("-vtop", "--vertical_top", nargs=1, default=['default'], help='Location of top of plot (vertical type) in mbar (pressure) or fractional height (height)')
parser.add_argument("-slay", "--split_layer", nargs=1, default=['no_split'], help='Split conserved quantities into weather and deep layers at this pressure')
parser.add_argument("-coord", "--coordinate_sys", nargs=1, default=['icoh'], help='For KE spectrum, use either icoh grid or llp grid')
parser.add_argument("-ladj", "--lmax_adjust", nargs=1, default=[0], help='For KE spectrum, icoh grid, adjust number of wave numbers to fit')
parser.add_argument("-slice", "--slice", nargs='+', default=['default'], help='Plot a long/lat slice or average over all values')
parser.add_argument("-mt", "--maketable", action='store_true', help='Print a table in text file of plot data')
parser.add_argument("-no-log", "--no_pressure_log", action='store_true', help='Switch off log coordinates in pressure (vertical)')
parser.add_argument("-llswap", "--latlonswap", action='store_true', help='Swap latitude and longitude axes (horizontal plots)')
parser.add_argument("-vc", "--vcoord", nargs=1, default=['pressure'], help='Vertical coordinate to use (pressure or height)')
parser.add_argument("-pgrid", "--pgrid_ref", nargs=1, default=['auto'], help='Reference file for pressure grid')
parser.add_argument("-clevs", "--clevels", nargs='+', default=[40], help='Set color contour levels')
args = parser.parse_args()

mph.make_plot(args)


last = time.time()
print(last - first)
