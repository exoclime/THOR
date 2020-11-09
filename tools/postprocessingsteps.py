#! /usr/bin/env python3
# -*- coding: utf-8 -*-

#---- code by Urs Schroffenegger ----------------------------------------------

"""Tool to rerun THOR on a folder and execute a small amount of postprocessing steps. 

Can rerun postprocessing steps configured by options:
* change resolution by using other opacities, spectrum and cloud files
* write out w0 and g0 values per band, layer, column
* write out f_dir spectrum per band, layer, column
* turn off Planck function
"""


import argparse

import pathlib

import datetime

import re
import subprocess
import sys
import shutil

def search_folder_for_pattern(folder, pattern):
    """Get a list of filenames and regex groups for filename matching pattern in folder."""
    files = folder.glob("*")

    return [(fn.name, reg.groups()) for fn, reg in
            map(lambda x : (x, re.compile(pattern).match(x.name)), files)
            if reg is not None] 

parser = argparse.ArgumentParser(description='execute THOR and alfrodull postprocessing step')

parser.add_argument("data_folder",
                    action='store',
                    default="./results/",
                    help="folder to look for data")

parser.add_argument('-f', '--postfix',
                    default=None,
                    action='store',
                    help="postfix to append to output dir")

parser.add_argument('-H', '--hires',
                    default=False,
                    action='store_true',
                    help="use hires spectrum")

parser.add_argument('--hiresfiles',
                    default="./Alfrodull/input/wasp43b/opac_sample_SI_r500.h5,./Alfrodull/input/stellar_spectrum_wasp43_r500.h5,./Alfrodull/input/clouds_enstatite_r500.h5",
                    action='store',
                    help="pathes to files used for hires spectrum, separated by commas")

parser.add_argument('-o', '--opacities',
                    default=False,
                    action='store_true',
                    help="output g0 and w0")

parser.add_argument('-b', '--beam',
                    default=False,
                    action='store_true',
                    help="output directional beam spectrum")

parser.add_argument('-p', '--noplanck',
                    default=False,
                    action='store_true',
                    help="set planck function to zero")

parser.add_argument('-n', '--numstep',
                    default=1,
                    type=int,
                    action='store',
                    help="number of postprocessing steps to run")

parser.add_argument('-c', '--numcol',
                    default=None,
                    type=int,
                    action='store',
                    help="number of parallel columns to run")

args = parser.parse_args()

base_folder = pathlib.Path(args.data_folder).resolve()

# search for last output file number
# get output write log
esp_write_log_list = search_folder_for_pattern(base_folder, r"^esp_write_log_.*\.txt$")

if len(esp_write_log_list) == 0: 
    print(f"Could not find output write log in folder {base_folder}")
    exit(-1)
elif len(esp_write_log_list) > 1: 
    print(f"Found more than one output write log in folder {base_folder}: {esp_write_log_list}")
    exit(-1)
else:
    print(f"Found output write log: {esp_write_log_list[0][0]}")
    esp_write_log = esp_write_log_list[0][0]

# get last written file
with (base_folder / esp_write_log).open("r") as f:
    lastline = f.readlines()[-1]
    (steps, f_index, f_name) = lastline.split()

print(f"Last output file {f_name} with index {f_index} at {steps} steps")

    
# search for last config file
all_config_copys = search_folder_for_pattern(base_folder, r"^config_copy\.(\d+)$")

# get last config copy
last_config_copy = sorted( all_config_copys, key=lambda x : x[1][0])[-1][0]

print(f"Found {len(all_config_copys)} configuration file copies, using last: {last_config_copy}")

# timestamp to use for this run
timestamp = datetime.datetime.now()
time_suffix = timestamp.strftime("%Y_%m_%d_%H_%M_%S")

# build list of replacements
regex_list = []
# replace number of steps
number_of_steps_to_add = int(args.numstep)

print(f"Running for {number_of_steps_to_add} step up to {int(steps) + number_of_steps_to_add}, storing output of every step")


regex_list.append((re.compile("^num_steps\s*=\s*(\d+)\s*", flags=re.MULTILINE), f"num_steps = {int(steps) + number_of_steps_to_add} "))
regex_list.append((re.compile("^n_out\s*=\s*(\d+)\s*", flags=re.MULTILINE), f"n_out = 1 "))

if args.opacities:
    print("Storing w0 and g0 per band")
    regex_list.append((re.compile("^Alf_store_w0_g0\s*=\s*((false)|(true))\s*", flags=re.MULTILINE), f"Alf_store_w0_g0 = true "))

if args.beam:
    print("Storing directional beam spectrum")
    regex_list.append((re.compile("^Alf_store_dir_spectrum\s*=\s*((false)|(true))\s*", flags=re.MULTILINE), f"Alf_store_dir_spectrum = true "))

if args.noplanck:
    print("Using null planck function")
    regex_list.append((re.compile("^Alf_null_planck_function\s*=\s*((false)|(true))\s*", flags=re.MULTILINE), f"Alf_null_planck_function = true "))
    
if args.hires:
    print("Running hires output")
    (opacity_filename, spectrum_filename, clouds_filename) = args.hiresfiles.split(",")
    print(f"Using opacity file: {opacity_filename}")
    print(f"Using spectrum file: {spectrum_filename}")
    print(f"Using clouds file: {clouds_filename}")
    regex_list.append((re.compile("^Alf_opacities_file\s*=\s*(.*)\s*", flags=re.MULTILINE), f"Alf_opacities_file = {opacity_filename}"))
    regex_list.append((re.compile("^Alf_stellar_spectrum\s*=\s*(.*)\s*", flags=re.MULTILINE), f"Alf_stellar_spectrum = {spectrum_filename} "))
    regex_list.append((re.compile("^Alf_cloudfile\s*=\s*(.*)\s*", flags=re.MULTILINE), f"Alf_cloudfile = {clouds_filename} "))

if args.numcol is not None:
    print("Setting number of parallel columns")
    regex_list.append((re.compile("^Alf_num_parallel_columns\s*=\s*(\d+)\s*", flags=re.MULTILINE), f"Alf_num_parallel_columns = {args.numcol} "))

if args.postfix is not None:
    postfix = args.postfix
else:
    postfix = time_suffix

with (base_folder / last_config_copy).open("r") as f:
    # read in the config
    config = f.read()
    # append signature
    config = f"# Configuration automatically updated by postprocessing script on {timestamp.isoformat(timespec='minutes')} \n\n" + config
    # perform the substitutions
    for rgx, sub in regex_list:
        print(f"Replacement configuration '{sub}'")
        (config, num_sub) = rgx.subn("\n#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#\n## REPLACED BY POSTPROCESSING SCRIPT\n"+sub + "\n#-#-#-#-#-#-#\n\n", config)
        if num_sub == 1:
            print("\tReplacement successfull")
        elif num_sub == 0:
            print("\tNo replacement found, appending")
            config = config + "\n\n#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#\n## INSERTED BY POSTPROCESSING SCRIPT\n"+  sub + "\n#-#-#-#-#-#-#\n\n"
        else:
            print("\tFound multiple entries for replacement")
            exit(-1)

    # print(config)

    
    post_processing_config = base_folder / f"config_postprocessing_{postfix}.thr"
    with post_processing_config.open("w") as w:
        w.write(config)
        


# Run thor on this config, wheeeeEEEEEEEEEEEE!
continue_from = base_folder / f_name
output_dir = base_folder / f"results_postprocessing_{postfix}"

esp_cmd = "./bin/esp"
esp_args = ["-c", f"{continue_from}", "-o", f"{output_dir}", f"{post_processing_config}"]
print(f"Running thor with:\n\t{' '.join([esp_cmd] + esp_args)}")
result = subprocess.run([esp_cmd] + esp_args,
                        stdout=sys.stdout,
                        stderr=sys.stderr,
                        universal_newlines=True)

if result.returncode != 0:
    print(f"Thor run returned bad exit code {result.returncode}")
    exit(-1)
    
    
# copy over grind and planet file
planet_files = search_folder_for_pattern(base_folder, r"^esp_output_planet_(.+).h5$")
if len(planet_files) == 1:
    # copy file over
    
    shutil.copy2(base_folder / planet_files[0][0], output_dir)
else:
    print(f"Found wrong number of planet file")
    for f, g in planet_files:
        print(f"\t[{f}]")
        
grid_files = search_folder_for_pattern(base_folder, r"^esp_output_grid_(.+).h5$")
if len(grid_files) == 1:
    # copy file over
    
    shutil.copy2(base_folder / grid_files[0][0], output_dir)
else:
    print(f"Found wrong number of grid file")
    for f, g in grid_files:
        print(f"\t[{f}]")
