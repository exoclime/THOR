
import configparser
import pathlib


# need command line arguments for those
base_output_dir = pathlib.Path('results_testing')

run_set_sel = 'fast'


# config sets
fast_set = [
    {'name': 'earth_hs',
     'base_ifile': 'ifile/earth_hstest.thr',
     'override': {'num_steps': '10'},
     'status': True,
     'compare_func': None,
     'compare_params': None},
    {'name': 'deephj',
     'base_ifile': 'ifile/deephj.thr',
     'override': {'num_steps': '10'},
     'status': True,
     'compare_func': None,
     'compare_params': None}
]

slow_set = [
    {'name': 'earth_hs',
     'base_ifile': 'ifile/earth_hstest.thr',
     'override': {'num_steps': '10000'},
     'status': True,
     'compare_func': None,
     'compare_params': None}
]


simulation_sets = {"slow": slow_set,
                   "fast": fast_set}


run_set = simulation_sets[run_set_sel]

if not base_output_dir.exists():
    base_output_dir.mkdir()
else:
    print("Output {} already exists, can't run".format(str(base_output_dir)))
    exit(-1)

test_results = {}

for config_set in run_set:
    print("Running {}".format(config_set['name']))

    config_parser = configparser.ConfigParser()
    f = open(config_set['base_ifile'])
    conf = "[config]\n" + f.read()
    config_parser.read_string(conf)

    # override configs
    for key, value in config_set['override'].items():
        config_parser['config'][key] = value

    output_dir = str(base_output_dir / config_set['name'])
    config_parser['config']['output'] = output_dir

    generated_config_name = base_output_dir / (config_set['name'] + ".thr")

    f = generated_config_name.open("w")

    for key, value in config_parser['config'].items():
        f.write("{} = {}\n".format(key, value))

    # run test

    # check output status

    # To be defined later, check output values if present
