
# THOR programming guide

## structure
### logic
#### dynamical core
#### physics
### Code structure
#### classes

## style
### definition
#### code and headers
     Each class should usually be in one header `*.h` and one code file `*.cpp`. Other than for test executables and main().
#### Alignement
##### indentation
      indent with tabs, 4 spaces for one tab.
##### line width
      Try to use 100 characters line width. Break long lines.
##### Braces 
      Open brace on same line as declaration, close it on its own line
##### Alignment
      - align similar vertical operators
      - align arguments
      - align comments
##### Naming
     CamelCase for class names, snake_case for functions and variables.
     

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
#### overview

#### binary comparison

#### tracing

## physical modules

