#---- Runs through a subset of simulations to test functionality of plot code -
#---- code by Russell Deitrick ------------------------------------------------

# to do: add a sync rot with pbl case
#        add hd 189 w/ and w/o full g
#        add W43 w/ HELIOS
#        clean up (delete unused output files or don't output at all)
#        put all in ifile/repo_benchmarks/
#        create script to run all

import subprocess as sub
import mjolnir_plot_helper as mph
import argparse
import os
from pathlib import Path
import time

#script designed to test mjolnir/regrid/hamarr functions and options
#things to add:
#   multi-panel plot
#   tidal earth with rotations
#   stream func fail with vc height
#   global diagnostics
#   specify options: slice, mt, clevs, vtop, stride (pgrid only? add to mjolnir?)

parser = argparse.ArgumentParser()
parser.add_argument('month_tag', nargs=1, help='Destination folder for figures (recommend using MonthYear like Jun2020)')
args = parser.parse_args()

# month_tag = 'Feb2022'  #change this to mandatory argument
month_tag = args.month_tag[0]
simulations_path = Path('testing/repo_benchmarks')
fig_destination_parent = Path('repo_benchmark_figures/'+month_tag)

sub_stdout = None

if not simulations_path.exists():
    raise IOError(simulations_path.__str__() +' does not exist!')
if not fig_destination_parent.exists():
    fig_destination_parent.mkdir()

# colors for console output
W = '\033[0m'  # white (normal)
R = '\033[1;31m'  # bold red
G = '\033[1;32m'  # bold green
O = '\033[33m'  # orange
B = '\033[34m'  # blue is hard to read on dark background
P = '\033[35m'  # purple
C = '\033[1;36m'  # bold cyan

def moveit(file,prefix,dest):
    #moves plots to new destination
    if not dest.exists():
        dest.mkdir()

    pname = file.split('/')[-1]
    sub.run(['cp '+file+' '+dest.__str__()+'/'+prefix+'_'+pname],shell=True)

def check_m_time(file,tol):
    mtime = os.path.getmtime(file)
    currenttime = time.time()
    if currenttime - mtime < tol:
        #was file modified in last tol seconds
        strout = file+' (over)written: '+C+'PASS'+W
    else:
        strout = file+' (over)written: '+R+'FAIL'+W
    return strout

tests = []
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# held suarez
sim_path = simulations_path / 'earth_hs'

### force overwrite the regrid process (check that files get written)
sub.run([f'pgrid -i 9 -l 54 {sim_path.__str__()} -w'],shell=True,
        stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_9_54_1.txt',60))

### check that we can make plots specifying pgrid file
sub.run([f'regrid -i 54 -l 54 {sim_path.__str__()} -w -pgrid pgrid_9_54_1.txt'],
        shell=True,stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/regrid_height_Earth_54.h5',60))
tests.append(
    check_m_time(f'{sim_path.__str__()}/pgrid_9_54_1/regrid_Earth_54.h5',60))

### mjolnir command line test
sub.run([f'mjolnir -i 54 -f {sim_path.__str__()} Tver uver Tlonver'],
        shell=True,stdout=sub_stdout)
tests.append(check_m_time(
        f'{sim_path.__str__()}/figures/temperature_p_ver_i54_l54_lat-90+00-90+00.pdf',60))

### figures to save for future comparison
fig_destination = fig_destination_parent / 'EarthHS_example'
prefix = 'earth_hs'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [9]
args.last_file = [54]
args.no_pressure_log = True
args.pview = ['Tver', 'uver', 'eddyKE', 'eddyMomMerid']
plots = mph.make_plot(args)
for p in plots:
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
            f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Earth RT sim (check Tsurface, insolation)
sim_path = simulations_path / 'earth_rt_dc_g5'

### regrid (check propagation of overwrite flag to pgrid and unmask flag)
sub.run([f'regrid -i 1 -l 1 {sim_path.__str__()} -w -unmask'],
        shell=True,stdout=sub_stdout)
#verify overwrite propagation
tests.append(check_m_time(
            f'{sim_path.__str__()}/pgrid_1_1_1.txt',60))
tests.append(check_m_time(
            f'{sim_path.__str__()}/pgrid_1_1_1/regrid_Earth_1.h5',60))

### check mjolnir plot helper thingie
fig_destination = fig_destination_parent / 'EarthRT_example'
prefix = 'earth_rt_unmask'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [1]
args.last_file = [1]
args.pview = ['all']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    if p == "'massf' plot type requires -vc height; plot not created":
        #this one should fail to produce a plot with these options
        tests.append(p+": "+C+"PASS"+W)
    else:
        moveit(p,prefix,fig_destination)
        tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

### remove regrid file, test mjolnir executes regrid, no unmask flag
os.remove(f'{sim_path.__str__()}/pgrid_1_1_1/regrid_Earth_1.h5')
prefix = 'earth_rt_mask'
args.pview = ['Tver', 'stream']
args.no_pressure_log = True
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
tests.append(check_m_time(
                f'{sim_path.__str__()}/pgrid_1_1_1/regrid_Earth_1.h5',60))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Wasp 43 b sim
sim_path = simulations_path / 'wasp43b_ex'

### force regrid
sub.run([f'regrid -i 59 -l 60 {sim_path.__str__()} -w'],shell=True,
    stdout=sub_stdout)

### test a few files to verify time averaging works
fig_destination = fig_destination_parent / 'WASP43b_example'
prefix = 'wasp43b'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [59]
args.last_file = [60]
args.vcoord = ['height']  #check height plotting/files
args.horizontal_lev = [1e3]
args.pview = ['uver','ulonver','Tulev','qheat','tracer']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

#check slice option works
args.slice = [-10,10]
args.pview = ['ulonver']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#hd189 w/ constant g and variable g
sim_path = simulations_path / 'hd189b_constg'

fig_destination = fig_destination_parent / 'HD189733b_constg'
prefix = 'hd189b_constg'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [30]
args.last_file = [30]
args.horizontal_lev = [10]
args.pview = ['uver', 'Tver', 'Tulev']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

sim_path = simulations_path / 'hd189b_fullg'

fig_destination = fig_destination_parent / 'HD189733b_fullg'
prefix = 'hd189b_fullg'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [30]
args.last_file = [30]
args.horizontal_lev = [10]
args.pview = ['uver', 'Tver', 'Tulev']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#sync rot pbl case
sim_path = simulations_path / 'sync_rot_pbl_test'

fig_destination = fig_destination_parent / 'sync_rot_pbl_test'
prefix = 'sync_rot_pbl'
sub.run([f'regrid -i 10 -l 60 {sim_path.__str__()} -rot -rot-ang 0 -90'],
        shell=True,stdout=sub_stdout)

args = mph.mjol_args(sim_path.__str__())
args.initial_file = [10]
args.last_file = [60]
args.horizontal_lev = [160]
args.no_pressure_log = True
args.pview = ['vver','wver','Kdiffver','Tulev','Tver','PTver','Fsens','Tsurf','stream']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

args.horizontal_lev = [20]
args.pview = ['Tulev']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#wasp43b w/ Alfrodull
sim_path = simulations_path / 'wasp43_ni_i2s'

fig_destination = fig_destination_parent / 'WASP43b_HELIOS'
prefix = 'wasp43b_helios'

args = mph.mjol_args(sim_path.__str__())
args.initial_file = [1]
args.last_file = [1]
args.horizontal_lev = [100]
args.pview = ['uver','Tver','Tulev','TSfluxprof','TSfdirprof','TSqheatprof','qheatprof']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(
                f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
print('\n'+G+"---"*20)
print("PLOT TEST FUNCTION REPORT")
print("---"*20+W+'\n')

for item in tests:
    print(item)
