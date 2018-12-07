
import configparser
import pathlib

import h5py
import numpy as np
import argparse

import asyncio


# argument parsing
parser = argparse.ArgumentParser(description='THOR tester tool')

parser.add_argument("-t",
                    "--testset",
                    action='store',
                    metavar='TESTSET',
                    help='test set to run',
                    default="fast")

parser.add_argument("-o",
                    "--output_dir",
                    action='store',
                    metavar='OUTPUTDIR',
                    help='directory for output of test files',
                    default="testing")

parser.add_argument("-d",
                    "--data_dir",
                    action='store',
                    metavar='DATADIR',
                    help='directory containing comparison data for test',
                    default="test_data")

args = parser.parse_args()

run_set_sel = args.testset

# need command line arguments for those
parent_test_dir = pathlib.Path(args.output_dir)
parent_data_dir = pathlib.Path(args.data_dir)

base_output_dir = parent_test_dir / run_set_sel
test_data_dir = parent_data_dir / run_set_sel

# colors for console output
W = '\033[0m'  # white (normal)
R = '\033[31m'  # red
G = '\033[32m'  # green
O = '\033[33m'  # orange
B = '\033[34m'  # blue is hard to read on dark background
P = '\033[35m'  # purple
C = '\033[1;36m'  # bold cyan

# HDF5 helper classes for comparison tests

def h5diff(ref_file, dat_file, epsilon=None):
    if not ref_file.exists():
        print("No ref file {}".format(ref_file))
        return False

    if not dat_file.exists():
        print("No dat file {}".format(dat_file))
        return False

    ref = h5py.File(str(ref_file))
    dat = h5py.File(str(dat_file))

    comp = True
    for datasetname in ref:
        if datasetname not in dat:
            print("No dataset {} in {}".format(datasetname, dat))
            comp = False
        else:
            # load data set
            dat_data = np.array(dat[datasetname])

            ref_data = np.array(ref[datasetname])
            if epsilon is None:
                if (dat_data != ref_data).any():
                    print("dataset {} mismatch in {}".format(datasetname, dat))
                    comp = False
            else:
                if (np.abs(dat_data - ref_data) > epsilon).any():
                    print("dataset {} mismatch in {} with epsilon {}".format(
                        datasetname, dat, epsilon))
                    comp = False

    return comp


def compare_h5files(name, output_dir, ref_dir, params=None):
    comparisons = params['comparisons']

    epsilon = None
    if 'epsilon' in params:
        epsilon = params['epsilon']

    equal = True

    for filename in comparisons:

        if not h5diff(ref_dir / name / filename, output_dir / name / filename, epsilon):
            equal = False
    # do some test, return value
    return equal

# Dummy test function for comparison


def testfunction(name, output_dir, ref_dir, params=None):
    # do some test, return value
    return True


######################################################################
# test sets definition
grid_and_startup_set = [
    # standard earth
    # test startup, grid file, planet file, and two first output files
    {'name': 'earth_hs_grid_4',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10',
                  'n_out': '10',
                  'glevel': '4',
                  'vlevel': '32'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Earth.h5',
                                        'esp_output_planet_Earth.h5',
                                        'esp_output_Earth_0.h5',
                                        'esp_output_Earth_1.h5']}},
    {'name': 'earth_hs_grid_5',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10',
                  'n_out': '10',
                  'glevel': '5',
                  'vlevel': '32'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Earth.h5',
                                        'esp_output_planet_Earth.h5',
                                        'esp_output_Earth_0.h5',
                                        'esp_output_Earth_1.h5']}},
    {'name': 'earth_hs_grid_6',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10',
                  'n_out': '10',
                  'glevel': '6',
                  'vlevel': '32'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Earth.h5',
                                        'esp_output_planet_Earth.h5',
                                        'esp_output_Earth_0.h5',
                                        'esp_output_Earth_1.h5']}},
]

# fast checks, just checking that it did not crash
fast_set = [
    # standard earth
    {'name': 'earth_hs',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10'},
     'status': 0,
     'compare_func': testfunction,
     'compare_params': {'param': 'novalue'}},

    # wrong config option, should fail
    {'name': 'shouldfail',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '00'},
     'status': 255,
     'compare_func': testfunction,
     'compare_params': {'param': 'novalue'}},

    # does not start from rest
    {'name': 'earth_hs_norest',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10',
                  'rest': 'false',
                  'initial': 'test_data/esp_initial.h5'},
     'status': 0,
     'compare_func': testfunction,
     'compare_params': {'param': 'novalue'}},


    # deepHJ
    {'name': 'deephj',
     'base_ifile': 'ifile/deephj.thr',
     'command_options': [],
     'override': {'num_steps': '10'},
     'status': 0,
     'compare_func': None,
     'compare_params': None},
    # Earth Sync
    {'name': 'earth_sync',
     'base_ifile': 'ifile/earth_sync.thr',
     'command_options': [],
     'override': {'num_steps': '100'},
     'status': 0,
     'compare_func': None,
     'compare_params': None},
    # ShallowHJ
    {'name': 'shallowhj',
     'base_ifile': 'ifile/shallowhj.thr',
     'command_options': [],
     'override': {'num_steps': '100'},
     'status': 0,
     'compare_func': None,
     'compare_params': None},
    # Planet of the Wasps
    {'name': 'wasp43b_ex',
     'base_ifile': 'ifile/wasp43b_ex.thr',
     'command_options': [],
     'override': {'num_steps': '100'},
     'status': 0,
     'compare_func': None,
     'compare_params': None},
]

# long tests
slow_set = [
    {'name': 'earth_hs',
     'base_ifile': 'ifile/earth_hstest.thr',
     'command_options': [],
     'override': {'num_steps': '10000',
                  'n_out': '2000' },
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Earth.h5',
                                       'esp_output_planet_Earth.h5',
                                       'esp_output_Earth_0.h5',
                                       'esp_output_Earth_1.h5',
                                       'esp_output_Earth_2.h5',
                                       'esp_output_Earth_3.h5'
                                       'esp_output_Earth_4.h5',
                                       'esp_output_Earth_5.h5' ]}},

     {'name': 'deephj',
     'base_ifile': 'ifile/deephj.thr',
     'command_options': [],
     'override': {'num_steps': '30000',
                  'n_out': '6000'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_deep.h5',
                                        'esp_output_planet_deep.h5',
                                        'esp_output_deep_0.h5',
                                        'esp_output_deep_1.h5',
                                        'esp_output_deep_2.h5',
                                        'esp_output_deep_3.h5'
                                        'esp_output_deep_4.h5',
                                        'esp_output_deep_5.h5' ]}},


     {'name': 'earth_sync',
     'base_ifile': 'ifile/earth_sync.thr',
     'command_options': [],
     'override': {'num_steps': '40000',
                  'n_out':'8000' },
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Earth.h5',
                                        'esp_output_planet_Earth.h5',
                                        'esp_output_Earth_0.h5',
                                        'esp_output_Earth_1.h5',
                                        'esp_output_Earth_2.h5',
                                        'esp_output_Earth_3.h5'
                                        'esp_output_Earth_4.h5',
                                        'esp_output_Earth_5.h5' ]}},

     {'name': 'shallowhj',
     'base_ifile': 'ifile/shallowhj.thr',
     'command_options': [],
     'override': {'num_steps': '40000',
                  'n_out': '8000'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_shallow.h5',
                                        'esp_output_planet_shallow.h5',
                                        'esp_output_shallow_0.h5',
                                        'esp_output_shallow_1.h5',
                                        'esp_output_shallow_2.h5',
                                        'esp_output_shallow_3.h5'
                                        'esp_output_shallow_4.h5',
                                        'esp_output_shallow_5.h5' ]}},

     {'name': 'wasp43b_ex',
     'base_ifile': 'ifile/wasp43b_ex.thr',
     'command_options': [],
     'override': {'num_steps': '10000',
                  'n_out': '2000'},
     'status': 0,
     'compare_func': compare_h5files,
     'compare_params': {'comparisons': ['esp_output_grid_Wasp43b.h5',
                                        'esp_output_planet_Wasp43b.h5',
                                        'esp_output_Wasp43b_0.h5',
                                        'esp_output_Wasp43b_1.h5',
                                        'esp_output_Wasp43b_2.h5',
                                        'esp_output_Wasp43b_3.h5'
                                        'esp_output_Wasp43b_4.h5',
                                        'esp_output_Wasp43b_5.h5' ]}},
]
######################################################################
# the simulation sets we can choose from
simulation_sets = {'slow': slow_set,
                   'fast': fast_set,
                   'grid': grid_and_startup_set}

run_set = simulation_sets[run_set_sel]

# make output directory
if not parent_test_dir.exists():
    parent_test_dir.mkdir()

if not base_output_dir.exists():
    base_output_dir.mkdir()
else:
    print("Output {} already exists, can't run".format(str(base_output_dir)))
    exit(-1)


# store results output for summary
test_results = {}


def log_result(name, result):
    print(name + ":\t" + result)
    if name in test_results:
        test_results[name].append(result)
    else:
        test_results[name] = [result]

# asynchronous function to run a process, capture its output and print out its output


async def run_subprocess(process_args):
    code = 'import datetime; print(datetime.datetime.now())'

    # Create the subprocess; redirect the standard output
    # into a pipe.
    proc = await asyncio.create_subprocess_exec(*process_args,
                                                stdout=asyncio.subprocess.PIPE,
                                                stderr=asyncio.subprocess.PIPE)

    # Read one line of output.
    result = ''
    while True:
        line = await proc.stdout.readline()
        if proc.stdout.at_eof():
            break
        if line != b'':
            print(" " + line.decode('ascii'), end='')
            result += line.decode('ascii')

    # Wait for the subprocess exit.
    await proc.wait()
    return result, proc.stderr, proc.returncode

# start event loop
loop = asyncio.get_event_loop()

for config_set in run_set:
    print(C+"Running {} in output directory {}".format(config_set['name'], str(base_output_dir))+W)

    config_parser = configparser.ConfigParser()
    config_parser.optionxform = lambda option: option
    f = open(config_set['base_ifile'])
    conf = "[config]\n" + f.read()
    config_parser.read_string(conf)

    # override configs
    for key, value in config_set['override'].items():
        config_parser['config'][key] = value

    output_dir = str(base_output_dir / config_set['name'])

    config_parser['config']['results_path'] = output_dir

    generated_config_name = base_output_dir / (config_set['name'] + ".thr")

    f = generated_config_name.open("w")

    for key, value in config_parser['config'].items():
        f.write("{} = {}\n".format(key, value))
    f.close()

    # run test
    command_options = config_set['command_options']

    print("starting bin/esp on {} with options {}".format(str(generated_config_name), command_options))
    stdout, stderr, returncode = loop.run_until_complete(run_subprocess(['bin/esp',
                                                                         str(generated_config_name)] + command_options
                                                                        ))

    # store output somewhere

    # check output status
    if returncode == config_set['status']:
        log_result(config_set['name'], G+"Finished running {} ended correctly".format(
            config_set['name'], returncode)+W)

        # check output data if we have a result evaluation function
        if config_set['compare_func'] is not None:
            compare_func = config_set['compare_func']
            compare_parameters = config_set['compare_params']
            compare_result = compare_func(
                config_set['name'], base_output_dir, test_data_dir, compare_parameters)
            if compare_result:
                log_result(config_set['name'], G+"data check passed"+W)
            else:
                log_result(config_set['name'], R+"data check failed"+W)

    else:
        log_result(config_set['name'], R+"Finished running {} failed with return code: ".format(
            config_set['name'], returncode) + W)
        log_result(config_set['name'], "return status for {}: {}".format(
            config_set['name'], returncode))
        log_result(config_set['name'], "stdout:\n {}".format(stdout))
        log_result(config_set['name'], "stderr:\n {}".format(stderr))


loop.close()
print(72*"*")
print(72*"*")
print("****     SUMMARY")
print(72*"*")
for name, result in test_results.items():
    for l in result:
        print(name + ":\t" + l)
