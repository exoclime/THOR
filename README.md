# THOR #

### Flexible Global Circulation Model to Explore Planetary Atmospheres

*THOR* is a GCM that solves the three-dimensional non-hydrostatic Euler equations on an icosahedral grid. *THOR* was designed to run on Graphics Processing Units (GPUs).

If you use this code please cite: [Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016](http://iopscience.iop.org/article/10.3847/0004-637X/829/2/115/meta)

###### Copyright (C) 2017 João Mendonça ######

### BUILD & RUN THOR

Main instructions to compile and run *THOR*. This version uses only a __single__ __GPU__.

Current code owner: Joao Mendonca: joao.mendonca@space.dtu.dk

Home website: http://software-oasis.com/

### INSTALL

*UBUNTU* *17.04*

1- Install cuda. 

```sh
   $ sudo apt install nvidia-cuda-toolkit
```
2- Downgrade g++ because cuda 8 conflicts with the latest g++.

```sh
   $ sudo apt install g++-5
```

3- Restart your pc.

4- Install hdf5 (https://support.hdfgroup.org/HDF5/release/obtainsrc.html). 
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

5- Go to THOR home directory.
   Open Makefile file and set "hdf_install" and "h5lib" paths correctly.

6- Set the value SM correctly, this depends on the GPU you are using. SM stands for Streaming Multiprocessor and the number indicates the features supported by the architecture. See https://developer.nvidia.com/cuda-gpus.
   Example: Tesla K20 -> 35. To get information on your GPU, type in the terminal: 

```sh
   $ nvidia-smi 
```
   
7- Create a config file in "/etc/ld.so.config.d" called for example "mylib.conf" with the following line:

   > /path_to_hdf5_libs
   
   Run: 
   
```sh
   $ sudo ldconfig 
```

8- Type: 

```sh
   $ make -j8
```

### RUN

*UBUNTU* *17.04*

1- Set planet's parameters in "src/initial/planet.cu".

2- Set model's parameter in "src/headers/define.h".

3- Run 

```sh
   $ ./bin/esp
```

4- Press enter and relax.


### Results

* Output is written in "results" folder.
* Very useful command lines to quickly explore the hdf5 files can be found in support.hdfgroup.org/HDFS/Tutor/cmdtools.html
  or type the command ">> man h5dump".
* You can find some Matlab routines to explore the results in "mjolnir" folder.
