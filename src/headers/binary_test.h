// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Description: binary correctness test of output, enabled with compile time switches
//
//   
//
// Method: [1] - Dumps output to binary file on a flag
//         [2] - Reads data from binary files on a flag and compare to
//               dumped output
//
// Known limitations: None.
//      
//
// Known issues: None.
//   
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// If you use this code please cite the following reference: 
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016  
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

#include "debug.h"


// precompiler macros to enable binary testing of simulation output
// BENCHMARKING enables the binary testing. if not set, those
// functions are empty and ignored by the compiler.
// BENCH_POINT_WRITE enables writing data out to reference files
// BENCH_POINT_COMPARE enables comparing current value with reference files
#ifdef BENCHMARKING
  #warning "Compiling with benchmarktest enabled"
  #define USE_BENCHMARK() binary_test & btester = binary_test::get_instance();
  #ifdef BENCH_POINT_WRITE // write reference mode
    #define BENCH_POINT(esp, iteration, name)  btester.output_reference(esp, iteration, name);
    #define BENCH_POINT_I(esp, iteration, name)  btester.output_reference(esp, std::to_string(iteration), name);
    #define BENCH_POINT_I_S(esp, iteration, subiteration, name)  btester.output_reference(esp, \
                                                                                          std::to_string(iteration) \
                                                                                          +"-" \
                                                                                          +std::to_string(subiteration), \
                                                                                          name);
    #define BENCH_POINT_GRID(grid)  btester.output_reference_grid(grid);
  #endif
  #ifdef BENCH_POINT_COMPARE // compare mode
    #define BENCH_POINT(esp, iteration, name)  btester.compare_to_reference(esp, iteration, name);
    #define BENCH_POINT_I(esp, iteration, name)  btester.compare_to_reference(esp, std::to_string(iteration), name);
    #define BENCH_POINT_I_S(esp, iteration, subiteration, name)  btester.compare_to_reference(esp, \
std::to_string(iteration)                                               \
                                                                                              +"-" \
                                                                                              +std::to_string(subiteration), \
                                                                                              name);
     #define BENCH_POINT_GRID(grid)  btester.compare_to_reference_grid(grid);
  #endif
#else // do nothing
#define USE_BENCHMARK()
  #define BENCH_POINT(esp, iteration, name)
  #define BENCH_POINT_I(esp, iteration, name)
  #define BENCH_POINT_I_S(esp, iteration, subiteration, name)
  #define BENCH_POINT_GRID(grid)
#endif

#ifdef BENCHMARKING
#include "esp.h"
#include "grid.h"
#include <memory>
#include "storage.h"


using std::string;

// singleton storing class for debug
class binary_test
{
public:
    // make a singleton, so that the object exists only once
    // use this to get a reference to the object
    static binary_test & get_instance();
    // no copy constructor and assignement operator
    binary_test(binary_test const&) = delete;
    void operator=(binary_test const&) = delete; 

    // esp object reference dump
    void output_reference(ESP & esp,
                          const string & iteration,
                          const string & ref_name);

    // esp comparison
    bool compare_to_reference(ESP & esp,
                              const string & iteration,
                              const string & ref_name);

    // grid object reference dump
    void output_reference_grid(Icogrid & grid);
    // grid object comparison
    bool compare_to_reference_grid(Icogrid & grid);

    void set_output(string base_name, string dir)
    {
        output_dir = dir;
        output_base_name = base_name;
    }
    
        
private:

    string output_dir;
    string output_base_name;
    
    
    // make constructor private, can only be instantiated through get_instance
    binary_test(string output_dir,
                string output_base_name);

    // helper function to compare two arrays for equality
    template<typename T>
    bool compare_arrays(int s1, T * d1,
                        int s2, T * d2,
                        string array = string(""),
                        bool print = false);

    // helper function to compare application array to saved array
    template<typename T>
    bool compare_to_saved_data(storage & s,
                               const string & name,
                               T * local_data,
                               const int & data_size);
  
};

// Compare binary table to saved table in storage output
template<typename T>
bool binary_test::compare_to_saved_data(storage & s,
                                        const string & name,
                                        T * local_data,
                                        const int & data_size) {
    std::unique_ptr<T[]> saved_data = nullptr;
    int size = 0;

//    cout << "Comparing " << name << " :\t";
    
    s.read_table(name, saved_data, size);

    
    
    bool b = compare_arrays(size, saved_data.get(),
                            data_size, local_data,
                            name, true);
//    cout << b << endl;
    return b;
    
}

// Binary comparison of two arrays d1 of size s1 and d2 of size s2
template<typename T>
bool binary_test::compare_arrays(int s1, T * d1,
                                 int s2, T * d2, string array, bool print)
{
    if (s1 != s2)
    {
        if (print)
            cout << array << ":\tdifferent sized arrays (" << s1
                 << ":" << s2 << endl;
        
        return false;
    }
    
    
    bool same = true;
    
    for (int i = 0; i < s1; i++)
        if (d1[i] != d2[i])
        {
            if (print)
                cout << array << "["<<i<<"]:\tdifferent value ("<<d1[i]<<":"<<d2[i]<<")"<<endl;
            
            same = false;
        }
    
    return same;
}
#endif // BENCHMARKING
