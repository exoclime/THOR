import numpy as np
import h5py
import asyncio
import pathlib
import configparser
import os
import scipy.interpolate as interp

# colors for console output
W = '\033[0m'  # white (normal)
R = '\033[31m'  # red
G = '\033[32m'  # green
O = '\033[33m'  # orange
B = '\033[34m'  # blue is hard to read on dark background
P = '\033[35m'  # purple
C = '\033[1;36m'  # bold cyan

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

def gen_h5_files(config_set):
    base_output_dir = pathlib.Path('ic_special')
    if not base_output_dir.exists():
        base_output_dir.mkdir()

    # start event loop
    loop = asyncio.get_event_loop()

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

    # check output status
    if returncode == config_set['status']:
        print('Generated h5 files in directory {}'.format(output_dir))
        os.system('mv {} {}'.format(output_dir+'/esp_output_'+config_parser['config']['simulation_ID']+'_0.h5',
                                    output_dir+'/esp_initial.h5'))
        os.system('mv {} {}'.format(output_dir+'/esp_output_planet_'+config_parser['config']['simulation_ID']+'.h5',
                                    output_dir+'/esp_initial_planet.h5'))
        os.system('rm {}'.format(output_dir+'/esp_*'+config_parser['config']['simulation_ID']+'*'))
        # add some grid etc settings for the editing step
        config_set['simulation_ID'] = config_parser['config']['simulation_ID']
        config_set['vlevel'] = config_parser['config']['vlevel']
        config_set['glevel'] = config_parser['config']['glevel']
        config_set['Top_altitude'] = config_parser['config']['Top_altitude']
    else:
        print('Something went wrong in THOR!')
    loop.close()

def edit_init_file(config_set):
    base_output_dir = pathlib.Path('ic_special')
    if not base_output_dir.exists():
        print("No *.h5 files have been generated in {}".format(str(base_output_dir)))
        exit(-1)

    output_dir = str(base_output_dir / config_set['name'])

    # here we load the esp_initial.h5 and edit the thermodynamic properties
    init_file = output_dir+'/esp_initial.h5'

    if os.path.exists(init_file):
        openh5 = h5py.File(init_file)
    else:
        print("Initial h5 file {} not found!".format(init_file))
        exit(-1)

    #--------vertical------------------------------
    nv = np.int(config_set['vlevel'])
    glev = np.int(config_set['glevel'])
    point_num = 2+10*2**(2*glev)
    dh = np.float(config_set['Top_altitude'])/nv
    final_height = np.arange(dh/2,np.float(config_set['Top_altitude']),dh)

    fvert = open(config_set['vertical_file'],'r+')
    vert_keys = fvert.readlines()[0].split()
    fvert.close()

    if not 'Height' in vert_keys:
        print("Error! Vertical data file {} must include 'Height'".format(config_set['vertical_file']))
        exit(-1)

    columns = np.loadtxt(config_set['vertical_file'],skiprows=1)
    vert_struct = {}
    for col in np.arange(len(vert_keys)):
        vert_struct[vert_keys[col]] = columns[:,col]

    vert_out = {}
    for key in vert_struct.keys():
        if key in openh5.keys():
            vert_out[key] = interp.pchip_interpolate(vert_struct['Height'],vert_struct[key],final_height)

    for i in np.arange(np.int(point_num)):
        for key in vert_out.keys():
            openh5[key][i*nv:i*nv+nv] = vert_out[key]

    openh5.close()
