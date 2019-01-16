#
# example of how to use thor_ic_gen to generate and edit initial conditions
# this sets up the initial conditions for earth_gwave_test.thr
#

import thor_ic_gen as ic

base_file = {'name': 'gen_ic_earth',
             'base_ifile': 'ifile/earth_gwave_test.thr',
             'command_options': [],
             'override': {'num_steps': '1',
                          'n_out': '10',
                          'glevel': '5',
                          'vlevel': '20',
                          'rest': 'true'},
             'status': 0,
             'vertical_file': 'ifile/TPprof_Earth_N1e-2.dat'}

#run to make thor generate the h5 files
ic.gen_h5_files(base_file)

#now edit the h5 files
ic.edit_init_file(base_file)
