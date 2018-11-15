
# THOR programming guide

## structure
### logic
The code uses a main loop, updating the state variables, going through these steps.
* dynamical core (THOR), updating all the atmospheric dynamics of the simulation
* physics (ProfX), 
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


you can get information on compilation error by compiling with `clang`
* error checking
```
$ make release -j8 COMP=clang++-6
```

This uses another compiler which repots a lot of issues. It also dumps a compilation database 
(in `compilation_commands.json`) usable by the clang tools (like clang-tidy) to check for errors 
and run sanitizers on the code. It can also be used by some editors to know how the files are compiled and help for navigation through code and code completion (ycmd for vim, irony-mode and rtags for emacs, ide-clangd for Atom).

### debugging tools
#### Overview
Debug tools are configured from `src/headers/debug.h`, by defining precompiler macros. They then execute some tests at check points that are written throughout the code. 


The check points are in the main loop, in the Thor and ProfX functions, checking values between different steps of the simulation. As we have loops, small steps and big steps, we mark them with loops within loops, that are given as first argument.

Check points look like this:

```
BENCH_POINT_I_S(current_step, rk, "Compute_Temperature_H_Pt_Geff", vector<string>({}), vector<string>({"temperature_d", "h_d", "hh_d", "pt_d", "pth_d", "gtil_d", "gtilh_d"}))
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
* The intermediate states go to `#define BENCHMARK_DUMP_REF_PATH` subdirectory
* compile
* run 
this saves reference files to the reference path. 

Then, apply the code changes you want to test, and enable dump comparison:
* enable debug tools, in `src/headers/debug.h`, uncomment `#define BENCHMARKING` 
* enable comparison of intermediate states. Uncomment `#define BENCHMARK_POINT_COMPARE`
* The intermediate states compare to the value of `#define BENCHMARK_DUMP_REF_PATH` subdirectory
* compile
* run 

This will print out data when the stored value isn't exactly equal to the computed value. 
It prints the bench-point name, the iteration numbers, the short name of arrays and 1 if the array matched, 0 if they didn't match, NA if it didn't find the array in the input file.


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



To have more verbose output, define `BENCH_PRINT_DEBUG`. This prints out all the comparisons, and not only when the comparison fail.

TODO: explain path thingy after implementation
TODO: explain statistics





#### tracing
The debug tools can help for debugging and tracing, by running checks and printing out results at each check points. The checks can be enabled in `debug.h` with  `BENCHMARKING` and the check flags.
* `BENCH_NAN_CHECK`: runs nan checks on the list of input arrays and prints out arrays that contain NaNs. 
* `BENCH_CHECK_LAST_CUDA_ERROR`: checks for cuda error codes and prints them out with the name of the 

You can add more checks by adding flags in `src/headers/debug.h` and the code in `src/devel/binary_test.cpp` in the function `binary_test::check_data()` with an appropriate `#ifdef` clause. 
   The function gets some text for the iteration count and the name of the steps, and a vector of tables to work on. It computes the definitions of tables (with pointers and sizes) to work on at the top.


## physical modules
New physical modules can be plugged into the dynamical core. 

## config file
### structure
### keys

## Running 
