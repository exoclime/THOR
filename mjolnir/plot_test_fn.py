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

simulations_path = Path('simulations/repo_benchmarks')
fig_destination = Path('simulations/repo_benchmarks/test_figures')

sub_stdout = None

if not simulations_path.exists():
    raise IOError(simulations_path.__str__() +' does not exist!')
if not fig_destination.exists():
    fig_destination.mkdir()

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

# held suarez
sim_path = simulations_path / 'earth_hs'
### force overwrite the regrid process (check that files get written)
sub.run([f'pgrid -i 9 -l 54 {sim_path.__str__()} -w'],shell=True,stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_9_54_1.txt',60))
sub.run([f'regrid -i 54 -l 54 {sim_path.__str__()} -w -pgrid pgrid_9_54_1.txt'],shell=True,stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/regrid_height_Earth_54.h5',60))
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_9_54_1/regrid_Earth_54.h5',60))
### mjolnir command line
sub.run([f'mjolnir -i 54 -f {sim_path.__str__()} Tver uver Tlonver'],shell=True,stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/figures/temperature_p_ver_i54_l54_lat-90+00-90+00.pdf',60))

# tidal earth to test rotations, etc

# Earth RT sim (check Tsurface, insolation)
sim_path = simulations_path / 'earth_rt_dc_g5'
### regrid (check propagation of overwrite flag to pgrid and unmask flag)
sub.run([f'regrid -i 60 -l 60 {sim_path.__str__()} -w -unmask'],shell=True,stdout=sub_stdout)
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_60_60_1.txt',60))  #verify overwrite propagation
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_60_60_1/regrid_Earth_60.h5',60))
### check mjolnir plot helper thingie
prefix = 'earth_rt_unmask'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [60]
args.last_file = [60]
args.pview = ['all']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    if p == "'massf' plot type requires -vc height; plot not created":
        #this one should fail to produce a plot with these options
        tests.append(p+": "+C+"PASS"+W)
    else:
        moveit(p,prefix,fig_destination)
        tests.append(check_m_time(f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

### remove regrid file, test mjolnir executes regrid, no unmask flag
os.remove(f'{sim_path.__str__()}/pgrid_60_60_1/regrid_Earth_60.h5')
prefix = 'earth_rt_mask'
args.pview = ['Tver', 'stream']
args.no_pressure_log = True
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
tests.append(check_m_time(f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
tests.append(check_m_time(f'{sim_path.__str__()}/pgrid_60_60_1/regrid_Earth_60.h5',60))

# Wasp 43 b sim (tracers??? someday!)
sim_path = simulations_path / 'wasp43b_example'
### force regrid
sub.run([f'regrid -i 59 -l 60 {sim_path.__str__()} -w'],shell=True,stdout=sub_stdout)
### test a few files to verify time averaging works
prefix = 'wasp43b'
args = mph.mjol_args(sim_path.__str__())
args.initial_file = [59]
args.last_file = [60]
args.vcoord = ['height']  #check height plotting/files
args.lev = [1e6]
args.pview = ['uver','ulonver','Tulev','qheat']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))
#check slice option works
args.slice = [-10,10]
args.pview = ['ulonver']
plots = mph.make_plot(args)
for p in plots:   # move to destination folder for easy viewing
    moveit(p,prefix,fig_destination)
    tests.append(check_m_time(f'{fig_destination.__str__()}/{prefix}_{p.split("/")[-1]}',60))

print('\n'+G+"---"*20)
print("PLOT TEST FUNCTION REPORT")
print("---"*20+W+'\n')

for item in tests:
    print(item)
