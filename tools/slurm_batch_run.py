"""
Queue multiple tasks for slurm in batch mode.
"""

import subprocess
import re
import argparse

import configparser
import pathlib

# Load default configuration
# create a file slurm.cfg in the directory you run this script from
# or in ~/.config/THOR-gcm/slurm.cfg.
# Following *.ini file format, key value pairs without quotes
# keys in config file will override keys from script
# ----------------
# [DEFAULTS]
# # working directory where slurm is run from
# working_dir = /path/to/thor/directory/
# # email address to send slurm report to
# user_email = thoruser@thorsim.org
# # where to log data in
# log_dir = /path/to/log/dir
# # slurm ressource request
# gpu_key = gpu:gtx1080ti:1
# # slurm partition
# partition = gpu

config_data = {'working_dir': "~/THOR/",             # path to THOR directory
               'user_email': "thoruser@thorsim.org",
               'log_dir': "~/",                      # path to log slurm output
               'gpu_key': 'gpu:1',                   # slurm argument for gpu selection  ( --gres )
               'partition': 'production'}           # slurm argument for partition ( -p )


def load_config():
    global config_data
    config_file = None
    # open default pathes
    cwd_config = pathlib.Path("slurm.cfg")
    # find in home directory, expanding '~'
    user_config = pathlib.Path("~/.config/THOR-gcm/slurm.cfg").expanduser()
    sys_config = pathlib.Path("/etc/THOR-gcm/slurm.cfg")

    configs_to_try = [cwd_config, user_config, sys_config]

    for path in configs_to_try:
        if path.exists():
            config_file = path
            break

    if config_file is not None:
        conf = configparser.ConfigParser()

        conf.read(config_file)

        for k, v in conf['DEFAULTS'].items():
            # example to convert some keys from text to boolean
            if (k == 'some_boolean_value'
                    or k == 'some_other_boolean_value'):
                config_data[k] = conf.getboolean('DEFAULTS', k)
            else:
                # read text keys
                config_data[k] = v


load_config()

print("config: ")
for k, v in config_data.items():
    print("\t{}: {}".format(k, v))
#######################################################################


parser = argparse.ArgumentParser()
parser.add_argument('input_file', metavar='nview', nargs='*', help='Thor config file to run')
parser.add_argument("-jn", "--job_name", nargs=1, default=['default'], help='Job name for slurm scheduler')
parser.add_argument("-n", "--num_jobs", nargs=1, default=[2], type=int, help='Number of sequential jobs to run (integer)')
parser.add_argument("-p", "--prof", action="store_true", default=False, help='Run profiler on job')
parser.add_argument("-o", "--output", type=str, default=None, help='Output dir name')
args = parser.parse_args()
initial_file = args.input_file[0]
if args.job_name[0] == 'default':
    job_name = args.input_file
else:
    job_name = args.job_name[0]
num_jobs = args.num_jobs[0]

if args.prof:
    profiling = f"{job_name}.prof"
else:
    profiling = None

if args.output is not None:
    output_arg = f"-o {args.output}"
else:
    output_arg = {}


def start_esp(args, esp_command, esp_args, initial_file, profiling=None):
    batch_id_re = re.compile("Submitted batch job (\d+)\n")
    print(f"Batch job args: {args}")
    with subprocess.Popen(args,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True) as proc:
        if profiling is not None:
            prof_cmd = f"nvprof --export-profile {profiling}"
        else:
            prof_cmd = ''
        proc.stdin.write("#!/bin/sh\n")
        srun_command = f"srun {prof_cmd} {esp_command} {esp_args} {initial_file}\n"
        proc.stdin.write(srun_command)
        proc.stdin.close()
        procid_string = proc.stdout.read()
        print(proc.stderr.read())
        m = batch_id_re.match(procid_string)
        if m is not None:
            batch_id = m.group(1)
            print(f"Submitted batch with ID: {batch_id}\nCommand: {srun_command}")
            return True, batch_id
        else:
            print(f"Error reading batch id from return: {procid_string}\nCommand: {srun_command}")

            return False, -1


# Configure batch
working_dir = str(pathlib.Path(config_data['working_dir']).expanduser())

# Acceptable  time  formats  include  "minutes", "minutes:seconds",
# "hours:minutes:seconds", "days-hours",
# "days-hours:minutes" and "days-hours:minutes:seconds".
time_limit = "0-23:00"
mail = config_data['user_email']
log_dir = pathlib.Path(config_data['log_dir']).expanduser()
if not log_dir.exists():
    print(f"Creating log output directory: {log_dir}")
    log_dir.mkdir(parents=True, exist_ok=True)

output_file = str(log_dir / f"slurm-esp-{job_name}-%j.out")  # %j for job index


args = ['sbatch',
        '-D', working_dir,
        '-J', job_name,
        '-n', str(1),
        '--gres', config_data['gpu_key'],
        '-p', config_data['partition'],
        '--time', time_limit,
        '--mail-type=ALL',
        '--mail-user=' + mail,
        '--output=' + output_file,
        '--signal=INT@60']
esp_command = "bin/esp"
esp_args = f'-b {output_arg}'


last_success = None
last_id = None
for i in range(num_jobs):
    if last_id is None:
        last_success, last_id = start_esp(args, esp_command, esp_args, initial_file, profiling=profiling)
    elif last_success:
        last_success, last_id = start_esp(args + ['--dependency=afterany:{}'.format(last_id)], esp_command, esp_args, initial_file, profiling=profiling)
    else:
        print("Error queuing last command")
        exit(-1)
