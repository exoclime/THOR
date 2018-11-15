# THOR #

### Flexible Global Circulation Model to Explore Planetary Atmospheres

*THOR* is a GCM that solves the three-dimensional non-hydrostatic Euler equations on an icosahedral grid. *THOR* was designed to run on Graphics Processing Units (GPUs).

If you use this code please cite: [Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016](http://iopscience.iop.org/article/10.3847/0004-637X/829/2/115/meta)

###### Copyright (C) 2017-2018 Exoclimes Simulation Platform ######

### Changes since version 1

* Addition of "tidally-locked Earth", "shallow hot jupiter", and "deep hot jupiter" __benchmark__ tests (see [Heng, Menou, & Phillips, 2011](https://academic.oup.com/mnras/article/413/4/2380/962712))

* Inclusion of __grey radiative transfer__ (see [Mendonca, J.M., Malik, M., Demory, B.-O., & Heng, K., AJ, 155, 150, 2018](http://iopscience.iop.org/article/10.3847/1538-3881/aaaebc/meta))

* Inclusion of top-of-atmosphere __Rayleigh drag__ ("sponge layer") (see [Mendonca, J.M., Tsai, S.-M., Malik, M., Grimm, S.L., & Heng, K.](http://adsabs.harvard.edu/abs/2018arXiv180800501M))

* Addition of __"conservation"__ routines, which calculate energy, entropy, mass, and angular momentum

* Inclusion of __dry convective adjustment__ scheme for sub-grid scale convection (used, though not detailed, in [Mendonca, J.M., Malik, M., Demory, B.-O., & Heng, K., AJ, 155, 150, 2018](http://iopscience.iop.org/article/10.3847/1538-3881/aaaebc/meta))

* Addition of __insolation__ calculation allowing for arbitrary orbit/rotation state (used with radiative transfer scheme)

* Compilation process has been completely overhauled to be more reliable and flexible

* Initial conditions are no longer hard coded, but are now set in user generated configuration files

* Command line options have been added to allow more flexibility when running the model, particularly for restarting canceled/finished simulations

* Modular structure for additional physics has been put in place and is used for the grey radiative transfer scheme 

* Numerous performance and debugging modes have been implemented in the code

* Previous MATLAB plotting routines have been adapted for Python. A number of new plotting options have been added to the Python code which were unavailable in the MATLAB code

* Output files now contain additional information about model settings and quantities related to the additions described above

### BUILD & RUN THOR

Main instructions to compile and run *THOR*. This version uses only a __single__ __GPU__.

Current code owners: Joao Mendonca: joao.mendonca@space.dtu.dk, Russell Deitrick: russell.deitrick@csh.unibe.ch, Urs Schroffenegger: urs.schroffenegger@csh.unibe.ch

### INSTALL

Tested on *UBUNTU* *17.04* *Debian unstable*

1- First, ensure that you have `git`, `make`, `gcc`, and `g++` installed. If you would like to use `cmake` to build THOR instead of `make`, ensure that that is installed as well. On Ubuntu, these can be installed like so
```sh
   $ sudo apt install git make
   $ sudo apt-get install gcc g++
```

2- Install CUDA. In Ubuntu, this can be done from the command line:

```sh
   $ sudo apt-get install nvidia-cuda-toolkit
```
Alternatively, you can find and download from the web: https://developer.nvidia.com/cuda-downloads

Note that you may have to manually add the paths to the nvidia compiler and libraries to your environment file. See "Chapter 7: Post-Installation Actions" in the installation guide here: https://developer.download.nvidia.com/compute/cuda/10.0/Prod/docs/sidebar/CUDA_Installation_Guide_Linux.pdf

3- Install HDF5, from your package manager if possible or by hand (see below)
```sh
   $ sudo apt-get install libhdf5-dev libhdf5-100  libhdf5-serial-dev libhdf5-cpp-100
```
For the python plotting scripts, you will need h5py. You can install it with your OS package manager or pip:
```sh
   $ pip3 install h5py
```

4- Use git to clone this repository to a location you can find 6 months from now.

### COMPILE THOR
#### Find your SM number

This depends on the GPU you are using. SM stands for Streaming Multiprocessor and the number indicates the features supported by the architecture. See https://developer.nvidia.com/cuda-gpus.
   Example: Tesla K20 -> 35. To get information on your GPU, type in the terminal:

```sh
   $ nvidia-smi
```

(cmake will try to guess that for you, if you compile with Makefile, you need to set this).

Depending on how you installed the CUDA-toolkit, you may also have to install some additional utilities to use `nvidia-smi`:
```sh
   $ sudo apt-get install nvidia-utils-390
```

#### Using makefile


Set your SM number in the makefile.
```
SM:=30 # Streaming Multiprocessor version
```
Or as argument on the command line.

Compile, for release build (`-j8` uses 8 parallel threads, to speed up compilation):

```sh
   $ make release -j8
```

for debug build:

```sh
   $ make debug -j8
```

* TRICK *
You can specify the SM number on the command line:
```sh
   $ make release -j8 SM=35
```

default build builds release and SM=30
```sh
   $ make -j8
```

To show commands echoed
```sh
   $ make VERBOSE=1
```

If if fails, check the makefile variables output at the beginning of the compilation, it shows the variables and the path detected for h5, which is a common cause of issue during compilation.

##### Define a local configuration Makefile.
Copy `Makefile.conf.template` to `Makefile.conf`. This defines a local makefile configuration that wont be added to git. You can define the `SM` number in there, so that you don't need to modify the main makefile that can be overwritten when pulling from git or add it to the command line each time.

#### Physics modules
You can use your own physics modules by setting the path to the physics module in the local makefile configuration file `Makefile.conf`, see [How to add your own physics modules](physics_modules.org).

#### Alternative: Using CMake
Thor also provides a CMake script to try to configure the SM and find the libraries for CUDA and HDF5 automatically.
Create a build directory and move into it:
```sh
   $ mkdir build
   $ cd build
```
Generate the makefile (don't forget the two points at the end of the command)

```sh
   $ cmake ..
```

Compile:
```sh
   $ make -j8
```

To recompile, run the make command again, no need to rerun cmake as long as you don't add files or change the build process.

* TRICKS *

To specify the SM architecture:
```sh
   $ cmake -DSM=30 ..
```

For more verbosity to debug makefile by showing commands:
```sh
   $ make VERBOSE=1
```

### INSTALL HDF5 from source

1- Install hdf5 (https://support.hdfgroup.org/HDF5/release/obtainsrc.html).
   Download the source code.
   Follow all the steps from the instructions, for example:

```sh
   $ cd <install_directory>

   $ mkdir build

   $ cd build

   $ sudo /path_to_HDF5_source/configure

   $ sudo make

   $ sudo make check

   $ sudo make install

   $ sudo apt install hdf5-helpers
```


2- Create a config file in "/etc/ld.so.config.d" called for example "mylib.conf" with the following line:

   > /path_to_hdf5_libs

   Run:

```sh
   $ sudo ldconfig
```

### RUN

*UBUNTU* *17.04*

1- Configure the initial conditions. Template of configuration is in "ifile/earth.thr".
Copy that file where you'd like as an initial condition file. e.g.: "init/myplanet.thr"

1- Set the planet's and model's parameters in "init/myplanet.thr".

2- Run

```sh
   $ ./bin/esp init/myplanet.thr
```

3- Press enter and go grab a coffee. Or lunch.

* command line arguments *
Positional argument: config filename (e.g. init/myplanet.thr)

Keyword argument:
```sh
 -g / --gpu_id <N>             GPU_ID to run on
 -o / --output_dir <PATH>      directory to write results to
 -i / --initial <PATH>         initial conditions HDF5 filename
 -N / --numsteps <N>           number of steps to run
 -w / --overwrite              Force overwrite of output file if they exist
 -c / --continue <PATH>        continue simulation from this output file
 -b / --batch                  Run as batch
```

Keyword arguments supersede config file arguments.
If initial conditions path is given on the command line, it starts from there instead of from rest and ignores the 'rest' setting in the config file.

* -g / --gpu_id

Uses the GPU configured by parameter

* -o / --output_dir
Writes results to this directory. It will also scan the output directory to check for already existing files and run, continue or restart depending on options.

* -N / --numsteps <N>
Number of steps of simulation to run.

* -i / --initial <PATH>
Instead of starting from rest, use <PATH> as initial conditions, using the provided model parameters. Checks consistency of models parameter with planet and grid definition used in initial file and starts from 0.



* -w / --overwrite

if output directory already contains files, the simulation does not run and outputs a warning. This forces the simulation to overwrite existing files.

* -c / --continue <PATH>
Continues a simulation from an output file. Like `--initial`, but continues at the simulation step and time from the input file. This provides the possibility to restart the simulation from a specific step (to contniue simulation, debug, or run with some changes in some parameters).

* -b / --batch
Run in batch mode in output directory. It checks output directory for result files:
 - if none exist, start a simulation from scratch.
 - if some exist, look for last valid written file, and continue from that point.
Useful to run simulation on a cluster with a time limit. When the application gets the INT or TERM signal, it writes down the last simulation step to disk.
Launching the simulation from a batch script with `-b` in the same folder starts the simulation or continues from the last save point point.

* exclusive options:
`--initial` and `--continue` are mutually exclusive.
`--batch` and  `--continue` are mutually exclusive.

#### Running benchmark tests
The benchmark tests will disable the aditional physics and enable forcing on the simulation. 
Set the `core_benchmark` key in the config file
Available values are:
1) HeldSuarez: Held-Suarez test for earth
2) HSShallowHotJupiter: Benchmark test for shallow hot Jupiter 
3) HSDeepHotJupiter: Benchmark test for deep hot Jupiter 
4) HSTidallyLockedEarth: Benchmark test for tidally locked Earth
5) NoBenchmark: No benchmark test - enables external physics module


### SLURM Batch script
#### Simple Batch script
Simple batch script launching SLURM on THOR, in `/home/thoruser/THOR`, with job name `earth`, configuration file `ifile/earth.thr`, on 1 task and 1 gpu, with a time limit of 2 hours, in file `esp.job`:

```sh
#!/bin/bash
#SBATCH -D /home/thoruser/THOR/
#SBATCH -J earth
#SBATCH -n 1 --gres gpu:1
#SBATCH --time 0-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=thoruser@thormail.com
#SBATCH --output="/home/thoruser/slurm-esp-%j.out"

srun bin/esp ifile/config.thr
```

Launch it in the job queue as
```
$ sbatch esp.job
```

#### Multiple batches with timeout

On a slurm queue with a time limit on jobs, you can run into the issue that the the queue kills your job after some time and you need to restart it to continue. For this, you need to queue several consecutive jobs with dependencies from one to the next.

Base script `simple_slurm.sh`, like the simple script, but starting in batch mode and sending interrupt 60 seconds before the end:



```sh
#!/bin/bash

#SBATCH -D /home/thoruser/THOR/
#SBATCH -J earth
#SBATCH -n 1 --gres gpu:1
#SBATCH --time 0-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=thoruser@thormail.com
#SBATCH --output="/home/thoruser/slurm-esp-%j.out"
#SBATCH --signal=INT@60

srun bin/esp ifile/config.thr -b
```

Batching script, `slurm_batch.sh` queuing a list of simple batches to restart.


```sh
#!/bin/bash

IDS=$(sbatch simple_slurm.sh)
ID=${IDS//[!0-9]/}
echo "ID $ID"

for t in {1..7..1}
do
  IDS=$(sbatch --dependency=afterany:$ID simple_slurm.sh)
  ID=${IDS//[!0-9]/}
  echo "ID $ID"
done
```


See `slurm_batch_run.py` for a Python script doing the same in one script.


### Results

* Output is written in "results" folder or path configured in config file or given as command line.
* Very useful command lines to quickly explore the hdf5 files can be found in support.hdfgroup.org/HDFS/Tutor/cmdtools.html
  or type the command ">> man h5dump".
* You can find some Matlab and Python routines to explore the results in "mjolnir" folder.

#### Python plotting

The mjolnir (`THOR/mjolnir/`) plotting scripts are written for Python 3 (compatibility with Python 2 still needs to be tested). Most dependencies are pretty standard for science: `numpy`, `matplotlib`, and `scipy`. Additionally, you'll need to install `h5py`:

```
$ pip3 install h5py
```

or

```
$ conda install h5py
```

`h5py` is simply a Python package that allows easy interaction with hdf5 files (THOR's standard output format).

`mjolnir.py` is set up as an executable for command line, but you'll need to add the path to your environment file. In bash, add the line

```
export PATH="$PATH:<path to thor>/mjolnir"
```

to your ~/.bashrc or ~/.bash_profile file, replacing `<path to thor>` with the actual path to the repository on your system. Probably not the smartest or most pythonic way of setting this up but one thing at a time please. Once that is done, the command to make a plot looks like

```
$ mjolnir <options> <type of plot>
```

For example,

```
$ mjolnir -i 0 -l 10 -f awesome_results -s coolest_planet Tver
```

where `-i`, `-l`, `-f`, and `-s` are options flags and `Tver` is a plot type. The available options flags are

```sh
 -i / --initial_file <N>       number of first output file to open
 -l / --last_file <N>          number of last output file to open
 -f / --file <string>          folder containing results
 -s / --simulation_ID <string> name of planet (used in naming of output files)
 -p / --pressure_lev <N>       pressure level to use in horizontal plots (mbar units)
 -pmin / --pressure_min <N>    pressure minimum for vertical plots (mbar units)
 -slay / --split_layer <N>     splits conservation data into "weather" and "deep" layers at this pressure (mbar units)
```

mjolnir averages the data over time for the entire range of files read in. So with `-i 0` and `-l 10`, files 0-10 will all be read in, and the plotted quantities will be averaged over all 11 snapshots in time. If you want to plot one instant in time, just set `-i` and `-l` to the same value. The averaging process can get quite long because the data is interpolated in many ways before averaging. Be careful if you are passing mjolnir more than ~50 output files.

There are three basic types of plot mjolnir can make (plus a few other special ones): vertical, horizontal, and profile. Vertical plots are averaged zonally and temporally, resulting in contours on a latitude vs pressure grid. Horizontal plots are averaged only temporally, and plotted at a given pressure level on a latitude vs longitude grid. Profile plots show quantities along single columns or average columns, as a function of pressure.

The `-p` option is used only by the horizontal plot types and is simply the desired pressure level to be viewed. The `-pmin` option is used only by the vertical plot types and is just the lowest pressure level to be plotted. For vertical plots, mjolnir will plot a dashed line representing the maximum pressure at the top of the model--thus data plotted above this line requires some extrapolation and should be viewed skeptically.

Current vertical plot types are

```sh
 Tver                           time- and zonally-averaged temperature
 uver                           time- and zonally-averaged zonal wind speed
 wver                           time- and zonally-averaged vertical wind speed
 PTver                          time- and zonally-averaged  potential temperature
 PVver                          time- and zonally-averaged  potential vorticity
 stream                         time- and zonally-averaged mass streaming function
```

Current horizontal plot types are

```sh
 Tulev                          time-averaged temperature and horizontal wind along a pressure surface
 ulev                           time-averaged zonal and meridional winds along a pressure surface   
 PVlev                          time-averaged potential vorticity along a pressure surface
 RVlev                          time-averaged relative vorticity along a pressure surface
 tracer                         time-averaged molecular abundances along a pressure surface
```

Current profile plot types are

```sh
  TP                            temperature-pressure profiles drawn from all over the grid
  wprof                         vertical wind vs pressure averaged over 4 quadrants (useful for synchronous rotation)   
```
