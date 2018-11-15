
# THOR programming guide

## structure
### logic
#### dynamical core
#### physics
### Code structure
#### classes

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
     
## Example
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
    The code can be formatted with [clang-format-8](https://clang.llvm.org/docs/ClangFormat.html). There is a clang-format-8 config file `.clang-format` in the root directory of the project that should apply. The tool can be integrated in vim, emacs and [Atom clang-format](https://atom.io/packages/clang-format).

You can run it from the command line running 

```
$ clang-format-8 -i <files>
```

This formats all the files passed as argument in place (with the `-i` option).

I use clang-format-8 to use the most recent options, if you don't have it, you might need to comment some options out of `.clang-format`.


you can get information on compilation error by compiling with `clang`
* error checking
```
$ make release -j8 COMP=clang++-8
```

This uses another compiler which repots a lot of issues. It also dumps a compilation database 
(in `compilation_commands.json`) usable by the clang tools (like clang-tidy) to check for errors 
and run sanitizers on the code.

### debugging tools
#### Overview
Debug tools are configured from `src/headers/debug.h`, by defining precompiler macros. They then execute some tests at check points that are written throughout the code. 


The check points are in the main loop, in the Thor and ProfX functions, checking values between different steps of the simulation. As we have loops, small steps and big steps, we mark them with loops within loops, that are given as first argument.

* `BENCH_POINT_I`: iteration number.
* `BENCH_POINT_I_S`: Iteration number and subiteration number.
* `BENCH_POINT_I_SS`: Iteration number, subiteration and subsubiteration number.

The other arguments are:
* Name of the bench point
* Vector of array names to track that are used as input to the next functions (currently unused)
* vector of array names to track that are used as output of the previous functions.

```
BENCH_POINT_I_S(current_step, rk, "Compute_Temperature_H_Pt_Geff", vector<string>({}), vector<string>({"temperature_d", "h_d", "hh_d", "pt_d", "pth_d", "gtil_d", "gtilh_d"}))
```

It then runs various checks depending on the other debug flags enabled.


#### binary comparison
The debug tools can dump the intermediate state of arrays from the simulation for each timesteps. And then in next runs, compare the computed values to the ones in the stored files. This helps to find if some code changes changed the results of some computations. 


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


You can add more arrays from the dynamical core to check by putting it in the list defined in `src/devel/binary_test.cpp`.

To have more verbose output, define `BENCH_PRINT_DEBUG`. This prints out all the comparisons, and not only when the comparison fail.

TODO: explain path thingy after implementation
TODO: explain statistics







#### tracing
The debug tools can help for debugging and tracing, by running checks and printing out results at each check points. The checks can be enabled in `debug.h` with  `BENCHMARKING` and the check flags.
* `BENCH_NAN_CHECK`: runs nan checks on the list of input arrays and prints out arrays that contain NaNs. 
* `BENCH_CHECK_LAST_CUDA_ERROR`: checks for cuda error codes and prints them out with the name of the 

You can add more checks by adding flags in `src/headers/debug.h` and the code in `src/devel/binary_test.cpp` in the function `binary_test::check_data()`.

## physical modules


## config file
### structure
### keys

### Running 
