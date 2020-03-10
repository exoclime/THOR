"""
Queue multiple tasks for slurm 
"""

import subprocess
import re

working_dir = "/home/urss/prog/THOR-dev/"
job_name = "wasp_alf_alf"
initial_file = "ifile/alf_wasp_alf.thr"
# Acceptable  time  formats  include  "minutes", "minutes:seconds",
# "hours:minutes:seconds", "days-hours",
# "days-hours:minutes" and "days-hours:minutes:seconds".
time_limit = "0-12:0"
mail = "urs.schroffenegger@csh.unibe.ch"

profiling=None
#profiling = f"{job_name}.prof"


output_file = f"/home/urss/slurm-esp-{job_name}-%j.out"
args = ['sbatch',
        '-D', working_dir,
        '-J', job_name,
        '-n', str(1),
        '--gres', 'gpu:1',
        '-p', 'production',
        '--time', time_limit,
        '--mail-type=ALL',
        '--mail-user='+mail,
        '--output='+output_file,
        '--signal=INT@60']
esp_command = "bin/esp_alf"
#esp_args = f'-b -N 20 -o {job_name}'
esp_args = f'-b -o {job_name}'
num_jobs = 1

def start_esp(args, esp_command, esp_args, initial_file, profiling=None):
    batch_id_re = re.compile("Submitted batch job (\d+)\n")

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
        proc.stdin.write(f"srun {prof_cmd} {esp_command} {esp_args} {initial_file}\n")
        proc.stdin.close()
        procid_string = proc.stdout.read()
        #print(procid_string)
        m = batch_id_re.match(procid_string) 
        if m is not None:
            batch_id = m.group(1)
            print("Submitted batch with ID: {}".format(batch_id))
            return True, batch_id
        else:
            print("Error reading batch id from return: {}".format(procid_string))
            
            return False, -1

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
