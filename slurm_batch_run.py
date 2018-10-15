"""                                                                                                                                                                              
Queue multiple tasks for slurm in batch mode. 
"""

import subprocess
import re


def start_esp(args, esp_command, esp_args, initial_file):
    batch_id_re = re.compile("Submitted batch job (\d+)\n")

    with subprocess.Popen(args,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True) as proc:
        proc.stdin.write("#!/bin/sh\n")
        proc.stdin.write("srun {} {} {}\n".format(esp_command,
                                                  esp_args,
                                                  initial_file))
        proc.stdin.close()
        procid_string = proc.stdout.read()
        m = batch_id_re.match(procid_string)
        if m is not None:
            batch_id = m.group(1)
            print("Submitted batch with ID: {}".format(batch_id))
            return True, batch_id
        else:
            print("Error reading batch id from return: {}".format(procid_string))

            return False, -1


# Configure batch
working_dir = "/path/to/thor/directory/"
job_name = "myjobname"
initial_file = "ifile/initialfile.thr"
# Acceptable  time  formats  include  "minutes", "minutes:seconds",
# "hours:minutes:seconds", "days-hours",
# "days-hours:minutes" and "days-hours:minutes:seconds".
time_limit = "0-0:10"
mail = "thoruser@thorsim.org"
output_file = "/path/to/slurm-output-%j.out"  # %j for job index
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
        '--signal=INT@60']  # send interrupt signal 60 seconds before end of job
esp_command = "bin/esp"
esp_args = '-b'

num_jobs = 10

# launch num_jobs slurm jobs with dependencies to the last job
last_success = None
last_id = None
for i in range(num_jobs):
    if last_id is None:
        last_success, last_id = start_esp(
            args, esp_command, esp_args, initial_file)
    elif last_success:
        last_success, last_id = start_esp(
            args + ['--dependency=afterany:{}'.format(last_id)], esp_command, esp_args, initial_file)
    else:
        print("Error queuing last command")
        exit(-1)
