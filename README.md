# THOR #

### Flexible Global Circulation Model to Explore Planetary Atmospheres

*THOR* is a GCM that solves the three-dimensional non-hydrostatic Euler equations on an icosahedral grid. *THOR* was designed to run on Graphics Processing Units (GPUs).

If you use this code please cite: [Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016](http://iopscience.iop.org/article/10.3847/0004-637X/829/2/115/meta)

###### Copyright (C) 2017-2018 Exoclimes Simulation Platform ######

### BUILD & RUN THOR

Main instructions to compile and run *THOR*. This version uses only a __single__ __GPU__.

Current code owners: Joao Mendonca: joao.mendonca@space.dtu.dk, Russell Deitrick: russel.deitrick@csh.unibe.ch, Urs Schroffenegger: urs.schroffenegger@csh.unibe.ch

### INSTALL

Tested on *UBUNTU* *17.04* *Debian unstable*

1- Install cuda. 

```sh
   $ sudo apt-get install nvidia-cuda-toolkit
```
2- Downgrade g++ because cuda 8/9 conflicts with the latest g++ and asks for g++-5

```sh
   $ sudo apt-get install g++-5
```
3- Install HDF5, from your package manager if possible or by hand (see below)
```sh
   $ sudo apt-get install libhdf5-dev libhdf5-100  libhdf5-serial-dev libhdf5-cpp-100 python-h5py
```
The python package is for analysis scripts.

### COMPILE THOR
#### Find your SM number

This depends on the GPU you are using. SM stands for Streaming Multiprocessor and the number indicates the features supported by the architecture. See https://developer.nvidia.com/cuda-gpus.
   Example: Tesla K20 -> 35. To get information on your GPU, type in the terminal: 

```sh
   $ nvidia-smi 
```

(cmake will try to guess that for you, if you compile with Makefile, you need to set this).

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
```

Keyword arguments supersede config file arguments. 
If initial conditions path is given on the command line, it starts from there instead of from rest and ignores the 'rest' setting in the config file.

### Results

* Output is written in "results" folder or path configured in config file or given as command line.
* Very useful command lines to quickly explore the hdf5 files can be found in support.hdfgroup.org/HDFS/Tutor/cmdtools.html
  or type the command ">> man h5dump".
* You can find some Matlab and Python routines to explore the results in "mjolnir" folder.
