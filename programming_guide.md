# THOR programming guide

## structure
### logic
The code uses a main loop, updating the state variables, going through these steps.
* dynamical core (THOR), updating all the atmospheric dynamics of the simulation
* physics (ProfX)
* conservation computation
* output


### Code structure
* Main function in `esp.cu`. Contains the argument reading, the setup, the initial condition calls and the main loop.
* grid: class to build hexagonal grid structure and operators on it. At startup, fills all the data arrays needed forcomputation on the grid.
* ESP, main class holding the state variables and the update function for the dynamical core and the physics. 
 - `src/initial/esp_initial.cu`: memory allocation and initial conditions loading and calculation.
 - `src/output/esp_output.cu`: file output functions to store the simulation state.
 - `src/thor/`: dynamical core step execution (`thor_driver.cu`) and kernel functions,
 - `src/headers/dyn/`: dynamical core headers and header only kernels
 - `src/profx/profx_driver.cu`: physics step execution. 
 - `src/headers/phy/`: physics headers and header only kernels.
* physics modules: 
 - `src/physics/` example physics modules, can be replace by user's own physics modules in their own folder outside of the dynamical core git tree.
* auxiliary classes
 - `src/input`: helpers for config file parsing and command arguments reading.
 - `src/utils/`: utilities, iteration timer helper class.
 - `src/files`: helpers for directories and logging.
 - `src/devel/`: helper for debugging and checking data consistency after changes.
 - `src/test/`: test programs for helper classes, code experiments. 
 
## style
* code and headers: Each class should usually be in one header `*.h` and one code file `*.cpp`. Other than for test executables and main().
* indentation: indent with tabs, 4 spaces for one tab.
* line width
Try to use 100 characters line width. Break long lines.
* Braces 
Open brace on same line as declaration, close it on its own line
* Alignment
- align similar vertical operators
- align arguments
- align comments
* Naming
CamelCase for class names, snake_case for functions and variables.
     
### Example
```
/* curly bracket on same line as definition */
// 4 spaces indentation

header: 

class MyObject{  
public: 
    MyObject(int blarp_):
    my_blarp(blarp_) 
    {

    };
    // Brackets
    int check_blarp(int blarp);

private:
    int my_blarp = 0;
    
}

class SubObject : 
    public MyObject {
public:

private: 
   int my_blurp;
}

cpp:
#include "MyObject.h"

int MyObject::CheckBlarp(int blarp){
    if (blarp == my_blarp){
       printf("We have blarp %d == %d\n",     // format string
               blarp,                          // input blarp
               my_blarp);                      // own blarp
    } else {
        printf("No blarp\n");
    } 
}


cpp:
int main(int argc, char* argv[]){

    MyObject my_object;
}

```

### Automatic formatting
The code can be formatted with [clang-format-6](https://clang.llvm.org/docs/ClangFormat.html). There is a clang-format-6 config file `.clang-format` in the root directory of the project that should apply. The tool can be integrated in vim, emacs and [Atom clang-format](https://atom.io/packages/clang-format).

You can run it from the command line running 

```
$ clang-format-6 -i <files>
```

This formats all the files passed as argument in place (with the `-i` option).

### Error checking and code analysis
you can get more information on the code by compiling with `clang`. `clang` is another compiler used on linux that now supports part of CUDA syntax, and is known for it's display of compilation messages. No performance benchmarks have been performed on it, it's only experimental to get more error checking on the code. Note that the clang output can be very verbose and very strict to conformance to standard.
* compiling with clang.
```
$ make release -j8 COMP=clang++-6
```

It also dumps a compilation database (in `compilation_commands.json`), which is usable by the clang tools (like clang-tidy) to do syntax analysis, modernising and run sanitizers on the code. It can also be used by some editors to know how the files are compiled and help for navigation through code and code completion (ycmd for vim, irony-mode and rtags for emacs, ide-clangd for Atom).

### debugging tools
#### Overview
Debug tools are configured from `src/headers/debug.h`, by defining precompiler macros. They then execute some tests at check points that are written throughout the code. 


The check points are in the main loop, in the Thor and ProfX functions, checking values between different steps of the simulation. As we have loops, small steps and big steps, we mark them with loops within loops, that are given as first argument.

Check points look like this:

```
BENCH_POINT_I_S(current_step, rk, "Compute_Temperature_H_Pt_Geff", (), ("temperature_d", "h_d", "hh_d", "pt_d", "pth_d", "gtil_d", "gtilh_d"))
```

The arguments are:
  - 1, 2 or 3 arguments (depending on using `BENCH_POINT_I`, `BENCH_POINT_I_S`, `BENCH_POINT_I_SS` describing the level in the update loop. (first iteration, sub-iteration and sub-sub-iteration)
   - String argument describing the step in update loop. (e.g. "RK2": second Runge Kutta step, "Vertical_Eq": vertical equilibrium computation.)
   - vector of string describing the arguments that are inputs to next step of simulation (currently unused, but can get stored to make comparisons between input and outputs)
   - vector of string describing the arguments that are outputs of the previous step of the simulation (before the call to the debug function)




It then runs various checks depending on the other debug flags enabled.


#### binary comparison
The debug tools can dump the intermediate state of arrays from the simulation for each timesteps. And then in next runs, compare the computed values to the ones in the stored files. This helps to check consistency when making code changes that shouldn't impact computation. 


Dump the reference files:
* enable debug tools, in `src/headers/debug.h`, uncomment `#define BENCHMARKING` 
* enable writing of intermediate states. Uncomment `#define BENCHMARK_POINT_WRITE`
* the intermediate states go to the output directory, in the `ref` subdirectory.
* the output files go now to the output directory, in the `write` subdirectory.
* compile
* run 
this saves reference files to the reference path. 

Then, apply the code changes you want to test, and enable dump comparison:
* enable debug tools, in `src/headers/debug.h`, uncomment `#define BENCHMARKING` 
* enable comparison of intermediate states. Uncomment `#define BENCHMARK_POINT_COMPARE`
* the intermediate states go to the output directory, in the `ref` subdirectory.
* the output files go now to the output directory, in the `compare` subdirectory.
* compile
* run 

This will print out data when the stored value isn't exactly equal to the computed value. 
It prints the bench-point name, the iteration numbers, the short name of arrays and `1` if the array matched, `0` if they didn't match, `NA` if it didn't find the array in the input file.


You can add more arrays from the dynamical core to check by putting it in the list defined in `src/devel/binary_test.cpp::build_definitions( ESP & esp, Icogrid & grid)`:
```
// {named index, { pointer, size, long name, short name, on device }}
 {"Rho_d",         { esp.Rho_d,         esp.nv*esp.point_num,   "Density", "rho", true}},
```
- named index: the index used in the code to find the info for that element
- pointer: the pointer to the data (on device or host)
- size: the size of the table
- name: the name to display 
- short name: the name to use in short debug summary
- on device: boolean telling if the pointer points to data on device or on host. If the data is on the device, it will copy the data to host or work on it on the device if needed.


##### Additional debug output
To have more verbose output, define `BENCH_PRINT_DEBUG`. This prints out all the comparisons, and not only when the comparison fail.

You can also print out statistics on errors by defining the `BENCH_COMPARE_PRINT_STATISTICS`. This will print statistics on tables with errors, (covering the failed values, when we talk about max and mean)
* absolute deviation: max and mean
* relative deviation: max and mean
* reference value: max and mean
* absolute value: max and mean


##### Fuzzy compare
By defining `BENCH_COMPARE_USE_EPSILON` and an epsilon value for `BENCH_COMPARE_EPSILON_VALUE`, it will try to compare the relative difference `abs(val - ref)/abs(ref)` to the epsilon value instead of an exact comparison.

Exact comparison is useful for bitwise comparison, but if algorithms are changed, the exact value will change but can stay close to the reference. For these tests, the fuzzy compare can be useful.




#### tracing
The debug tools can help for debugging and tracing, by running checks and printing out results at each check points. The checks can be enabled in `debug.h` with  `BENCHMARKING` and the check flags.
* `BENCH_NAN_CHECK`: runs nan checks on the list of input arrays and prints out arrays that contain NaNs. 
* `BENCH_CHECK_LAST_CUDA_ERROR`: checks for cuda error codes and prints them out with the name of the 

You can add more checks by adding flags in `src/headers/debug.h` and the code in `src/devel/binary_test.cpp` in the function `binary_test::check_data()` with an appropriate `#ifdef` clause. 
   The function gets some text for the iteration count and the name of the steps, and a vector of tables to work on. It computes the definitions of tables (with pointers and sizes) to work on at the top.

## Programming tools
### Reduction Add
Some operations need to perform an addition over all elements of an array on the device. To execute this, instead of running some AtomicAdd that is suboptimal and not binary reproducible, you can run a kernel that parallelises the operation.

Use the reduction add operation defined in `src/headers/reduction_add.h`. To run a sum on an array on the device, use:

```
double gpu_sum_on_device<BLOCKSIZE>(double *data_d, long length);
```

Arguments:
* `BLOCKSIZE`: template argument for size of blocks to work on. Must be a power of two, 1024 is a good value. It will sum over sub-arrays of this size on the device and then sum up the results on the host, so bigger is better.
* `data_d`: pointer to memory on the device to operate on. This memory has to be aligned to memory boundaries. Some buffers needed padding to get to the correct alignement.
* `length`: size of the array to operate on. Doesn't need to be a power of 2, the function will do the padding. Can be smaller than `BLOCKSIZE`.


### 3D Vector operations 

To work on vector values, CUDA provides some basic N-dimensional data types. THOR can use `double3` for 3D vectors. It defines operators in `src/headers/vector_operations.h` so that this datatype can be used in amore convenient way. 

* standard math operators `+`,`-`, adding two vectors or a scalar to all elements of a vector. (also as inplace operators, `+=`, `-=`). 
* multiplication by a scalar `*` and division `/` by a scalar. (in place as `*=` and `/=`)
* dot product `dot(v1,v2)`
* length `length(v)`
* normalisation `normalize(v)`
* cross product `cross(v1, v2)`

Usage examples in `src/grid/grid.cu`.

Some versions of CUDA provide their own version of these operators, but not available on all platforms tested.

## physical modules
New physical behaviour can be implemented outside of the main dynamical core and plugged in as modules.

Two examples are available:
* `src/physics/managers/empty/` is an empty code structure to implement user modules
* `src/physics/managers/multi/` is an example of a physics modules using multiple physics behaviour, implemented in subclasses, that provide radiative transfer and chemical dynamics, with the implementation of the physics in itself  in `src/physics/modules`.

#### Interface
The physics modules themselves must provide the function signatures as defined in `src/headers/phy_modules.h`. 

They should mange their own state variables and functions, and fill in the hooks to be executed from the main loop.

- `phy_modules_get_name`: called at initialisation to store name of module used in output file, for reference.
- `phy_modules_print_config`: called at initialisation for CLI reporting of configuration.
- `phy_modules_generate_config`: called with `config_reader` to add configurations keys that should be read from input file.
- `phy_modules_init_mem`: called at initialisation for module to allocate it's memory and cuda memory, receivees a `device_RK_array_manager` to register the arrays that should be updated by the RK kernels during the dynamical core Runge-Kutta step.
- `phy_modules_free_mem`: called at the end of the application to free memory.
- `phy_modules_init_data`: called at initialisation to initialise the state variables. Receives main state variables and a pointer to astorage helper object. If storage object is null, starting from rest, initialise the module with default. If storage pointer is non null, wraps the start up file used to run thor, it can read it's own state in that file if it has been save to and restart from there.
- `phy_modules_dyn_core_loop_init`: called before the dynamical core step. Usually used to initialise data for the step or swap data from step initial state and step final state arrays.
- `phy_modules_dyn_core_loop_slow_modes`: called during slow step of dynamical core integration. 
- `phy_modules_dyn_core_loop_fast_modes`: called during fast step of dynamics core integrations.
- The arrays registered in `phy_modules_init_mem` are advanced in `UpdateRK` and `UpdateRK2` through aRunge-Kutta scheme.
- `phy_modules_dyn_core_loop_end`: end of dynamical core loop. Used to swap initial/final state arrays before physical modules step.

- `phy_modules_phy_loop`: physics integration scheme.
- `phy modules_store`: called at end of N step to store data from integration.
- `phy_modules_store_init`: called at initalisation to save parameters of physics module.

#### Compilation of module 
The main makefile calls the makefile in the physics module directory. It passes it its own variables to help for compilation. 

The physics module should create a static library called `libphy_modules.a` in its root directory that will be used to link to in the main program. 

#### Integration of module in THOR.
To integrate your physical module to THOR, configure it's relative path in `Makefile.conf`, by setting the `MODULES_SRC` variable.

The main makefile will then run the makefile in that directory and statically link to the `libphy_modules.a` library found at that path.

## Mesh structure and indexing



## config file
### structure
### keys

## Running 
