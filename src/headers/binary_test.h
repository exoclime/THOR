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

#include "../headers/debug.h"
#include "../headers/esp.h"
#include <memory>
#include "storage.h"

#ifdef BENCHMARKING
  #warning "Compiling with benchmarktest enabled"
#define USE_BENCHMARK(esp) binary_test & btester = binary_test::get_instance(esp);
  #ifdef WRITE_BINARY // write reference mode
    #define BENCH_POINT(iteration, name)  btester.output_reference(iteration, name);
  #else // compare mode
    #define BENCH_POINT(iteration, name) btester.compare_to_reference(iteration, name);
  #endif 
#else // do nothing
  #define USE_BENCHMARK(esp)
  #define BENCH_POINT(iteration, name)
#endif




using std::string;

enum class binary_test_mode {grid, data};

// singleton storing class for debug
class binary_test
{
public:
    static binary_test & get_instance(ESP & esp);
    binary_test(binary_test const&) = delete;
    void operator=(binary_test const&) = delete; 

    void output_reference(const int & iteration,
                          const string & ref_name,
                          const binary_test_mode & mode = binary_test_mode::data);

    bool compare_to_reference(const int & iteration,
                              const string & ref_name,
                              const binary_test_mode & mode = binary_test_mode::data);
    
private:
    ESP & esp;

    string output_dir;
    string output_base_name;
    
    
    // make constructor private, can only be instantiated through get_instance
    binary_test(ESP & esp_,
                string output_dir,
                string output_base_name);

    // helper function to compare two arrays for equality
    template<typename T>
    bool compare_arrays(int s1, T * d1,
                        int s2, T * d2);

    // helper function to compare application array to saved array
    template<typename T>
    bool compare_to_saved_data(storage & s,
                               const string & name,
                               T * local_data,
                               const int & data_size);
    
   
    
//    binary_test(binary_test const&);              // Don't Implement
//    void operator=(binary_test const&); // Don't implement
  
};

template<typename T>
bool binary_test::compare_to_saved_data(storage & s,
                                        const string & name,
                                        T * local_data,
                                        const int & data_size) {
    std::unique_ptr<T[]> saved_data = nullptr;
    int size = 0;
    
    s.read_table(name, saved_data, size);
    
    return compare_arrays(size, saved_data.get(),
                          data_size, local_data);
}
                          
template<typename T>
bool binary_test::compare_arrays(int s1, T * d1,
                                 int s2, T * d2)
{
    if (s1 != s2)
        return false;
    
    bool same = true;
    
    for (int i = 0; i < s1; i++)
        if (d1[i] != d2[i])
            same = false;
    return same;
}
