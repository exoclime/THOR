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
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <vector>

#include "debug.h"
#include "directories.h"
#include <algorithm>
#include <iomanip>
#include <memory>
// precompiler macros to enable binary testing of simulation output
// BENCHMARKING enables the binary testing. if not set, those
// functions are empty and ignored by the compiler.
// BENCH_POINT_WRITE enables writing data out to reference files
// BENCH_POINT_COMPARE enables comparing current value with reference files

#ifdef BENCHMARKING
#    define BRACED_INIT_LIST(...) \
        { __VA_ARGS__ }
#    define USE_BENCHMARK() binary_test& btester = binary_test::get_instance();
#    define SET_BENCHMARK_PATH(path) \
        binary_test& btester = binary_test::get_instance().set_output("bindata_", path);
#    define INIT_BENCHMARK(esp, grid, path)                                           \
        binary_test::get_instance().append_definitions(build_definitions(esp, grid)); \
        binary_test::get_instance().set_output("bindata_", path);
#    define BENCH_POINT(iteration, name, in, out)                          \
        btester.check_data(iteration,                                      \
                           name,                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),  \
                           std::vector<std::string>(BRACED_INIT_LIST out), \
                           false);
#    define BENCH_POINT_I(iteration, name, in, out)                        \
        btester.check_data(std::to_string(iteration),                      \
                           name,                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),  \
                           std::vector<std::string>(BRACED_INIT_LIST out), \
                           false);
#    define BENCH_POINT_I_S(iteration, subiteration, name, in, out)                        \
        btester.check_data(std::to_string(iteration) + "-" + std::to_string(subiteration), \
                           name,                                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),                  \
                           std::vector<std::string>(BRACED_INIT_LIST out),                 \
                           false);

#    define BENCH_POINT_I_SS(iteration, subiteration, subsubiteration, name, in, out)           \
        btester.check_data(std::to_string(iteration) + "-" + std::to_string(subiteration) + "-" \
                               + std::to_string(subsubiteration),                               \
                           name,                                                                \
                           std::vector<std::string>(BRACED_INIT_LIST in),                       \
                           std::vector<std::string>(BRACED_INIT_LIST out),                      \
                           false);

#    define BENCH_POINT_PHY(iteration, name, in, out)                      \
        btester.check_data(iteration,                                      \
                           name,                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),  \
                           std::vector<std::string>(BRACED_INIT_LIST out), \
                           true);
#    define BENCH_POINT_I_PHY(iteration, name, in, out)                    \
        btester.check_data(std::to_string(iteration),                      \
                           name,                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),  \
                           std::vector<std::string>(BRACED_INIT_LIST out), \
                           true);
#    define BENCH_POINT_I_S_PHY(iteration, subiteration, name, in, out)                    \
        btester.check_data(std::to_string(iteration) + "-" + std::to_string(subiteration), \
                           name,                                                           \
                           std::vector<std::string>(BRACED_INIT_LIST in),                  \
                           std::vector<std::string>(BRACED_INIT_LIST out),                 \
                           true);
#    define BENCH_POINT_I_SS_PHY(iteration, subiteration, subsubiteration, name, in, out)       \
        btester.check_data(std::to_string(iteration) + "-" + std::to_string(subiteration) + "-" \
                               + std::to_string(subsubiteration),                               \
                           name,                                                                \
                           std::vector<std::string>(BRACED_INIT_LIST in),                       \
                           std::vector<std::string>(BRACED_INIT_LIST out),                      \
                           true);


#else // do nothing
#    define USE_BENCHMARK()
#    define SET_BENCHMARK_PATH(path)
#    define INIT_BENCHMARK(esp, grid, path)
#    define BENCH_POINT(iteration, name, in, out)
#    define BENCH_POINT_I(iteration, name, in, out)
#    define BENCH_POINT_I_S(iteration, subiteration, name, in, out)
#    define BENCH_POINT_I_SS(iteration, subiteration, subsubiteration, name, in, out)
// also trace data from physics module
#    define BENCH_POINT_PHY(iteration, name, in, out)
#    define BENCH_POINT_I_PHY(iteration, name, in, out)
#    define BENCH_POINT_I_S_PHY(iteration, subiteration, name, in, out)
#    define BENCH_POINT_I_SS_PHY(iteration, subiteration, subsubiteration, name, in, out)
#endif // BENCHMARKING

#ifdef BENCHMARKING
#    include "esp.h"
#    include "grid.h"
#    include "storage.h"
#    include <functional>
#    include <map>
#    include <memory>
#    include <vector>

using std::map;
using std::string;
using std::vector;


struct output_def {
    double*&                             data;
    int                                  size;
    string                               name;
    string                               short_name;
    bool                                 device_ptr;
    std::function<std::string(int, int)> fun;
};

struct compare_statistics {
    int    num_values;
    int    num_failures;
    int    first_failure_idx;
    double max_abs_delta;
    double mean_abs_delta;
    double max_rel_delta;
    double mean_rel_delta;
    double mean_val;
    double max_val;
    double mean_ref;
    double max_ref;
};


map<string, output_def> build_definitions(ESP& esp, Icogrid& grd);

// singleton storing class for debug
class binary_test
{
public:
    // make a singleton, so that the object exists only once
    // use this to get a reference to the object
    static binary_test& get_instance();
    // no copy constructor and assignement operator
    binary_test(binary_test const&) = delete;
    void operator=(binary_test const&) = delete;

    // esp reference dump
    void output_reference(const string&             iteration,
                          const string&             ref_name,
                          const vector<output_def>& output_reference);

    // comparison
    bool compare_to_reference(const string&             iteration,
                              const string&             ref_name,
                              const vector<output_def>& output_reference);


    void check_data(const string&         iteration,
                    const string&         ref_name,
                    const vector<string>& input_vars,
                    const vector<string>& output_vars,
                    const bool&           trace_phy_module);

    void set_output(string base_name, string dir) {
        output_dir       = (path(dir) / string("ref")).to_string();
        output_crash     = (path(dir) / string("crash")).to_string();
        output_base_name = base_name;
    }


    void append_definitions(const std::map<string, output_def>& defs);

    // make constructor private, can only be instantiated through get_instance
    ~binary_test();

    void register_phy_modules_variables(const std::map<string, output_def>& defs,
                                        const vector<string>&               phy_modules_data_in,
                                        const vector<string>&               phy_modules_data_out);

private:
    vector<string> phy_modules_data_in;
    vector<string> phy_modules_data_out;

    vector<string> current_input_vars;


    string output_dir;
    string output_crash;
    string output_base_name;


    // make constructor private, can only be instantiated through get_instance
    binary_test(string output_dir, string output_base_name);

    // helper function to compare two arrays for equality
    template<typename T>
    bool compare_arrays(int                 s1,
                        T*                  d1,
                        int                 s2,
                        T*                  d2,
                        compare_statistics& stats,
                        string              array = string(""),
                        bool                print = false);

    // helper function to compare application array to saved array
    template<typename T>
    bool compare_to_saved_data(storage&            s,
                               const string&       name,
                               T*                  local_data,
                               const int&          data_size,
                               compare_statistics& stats);

    bool check_nan(const string&             iteration,
                   const string&             ref_name,
                   const vector<output_def>& data_output);


    bool                         output_defined = false;
    std::map<string, output_def> output_definitions;
    std::unique_ptr<double[]>    mem_buf = nullptr;

    bool* nan_check_d;
};

// Compare binary table to saved table in storage output
template<typename T>
bool binary_test::compare_to_saved_data(storage&            s,
                                        const string&       name,
                                        T*                  local_data,
                                        const int&          data_size,
                                        compare_statistics& stats) {
    std::unique_ptr<T[]> saved_data = nullptr;
    int                  size       = 0;
    // print out piece of the table for visual inspection
    bool print_details = false;

    //    cout << "Comparing " << name << " :\t";

    s.read_table(name, saved_data, size);


    bool b =
        compare_arrays(size, saved_data.get(), data_size, local_data, stats, name, print_details);
    //    cout << b << endl;
    return b;
}

// Binary comparison of two arrays d1 of size s1 and d2 of size s2
template<typename T>
bool binary_test::compare_arrays(int                 s1,
                                 T*                  d1,
                                 int                 s2,
                                 T*                  d2,
                                 compare_statistics& stats,
                                 string              array,
                                 bool                print) {
#    ifdef BENCH_COMPARE_PRINT_STATISTICS
    stats.num_values        = 0;
    stats.num_failures      = 0;
    stats.first_failure_idx = 0;
    stats.max_abs_delta     = 0.0;
    stats.mean_abs_delta    = 0.0;
    stats.max_rel_delta     = 0.0;
    stats.mean_rel_delta    = 0.0;
    stats.mean_val          = 0.0;
    stats.max_val           = 0.0;
    stats.mean_ref          = 0.0;
    stats.max_ref           = 0.0;
#    endif // BENCH_COMPARE_PRINT_STATISTICS

    if (s1 != s2) {
        if (print)
            log::printf(":\tdifferent sized arrays (%d:%d)\n", array.c_str(), s1, s2);

        return false;
    }


    bool same = true;
#    ifdef BENCH_COMPARE_PRINT_STATISTICS
    bool first_failure = true;
#    endif // BENCH_COMPARE_PRINT_STATISTICS

    for (int i = 0; i < s1; i++) {

        stats.num_values += 1;

        //>
        //double mx = (abs(d1[i]) > abs(d2[i]))?abs(d1[i]):abs(d2[i]);
#    ifdef BENCH_COMPARE_USE_EPSILON
        if (abs((d1[i] - d2[i])) > std::abs(d1[i]) * BENCH_COMPARE_EPSILON_VALUE) {
#    else  // compare absolute value
        if (d1[i] != d2[i]) {
#    endif // BENCH_COMPARE_USE_EPSILON

            //            if (print && i < 10)
            if (print)
                log::printf(
                    "%s [ %d ]:\tdifferent value (%.20e:%.20e)\n", array.c_str(), i, d1[i], d2[i]);

#    ifdef BENCH_COMPARE_PRINT_STATISTICS
            stats.num_failures += 1;

            if (first_failure) {
                stats.first_failure_idx = i;
                first_failure           = false;
            }


            T      delta     = std::abs(d1[i] - d2[i]);
            double delta_rel = delta / std::abs(d1[i]);

            stats.mean_val += std::abs(d2[i]);
            stats.max_val = std::max(d2[i], stats.max_val);

            stats.mean_ref += std::abs(d1[i]);
            stats.max_ref = std::max(d1[i], stats.max_ref);

            stats.mean_abs_delta += delta;
            stats.mean_rel_delta += delta_rel;
            stats.max_abs_delta = std::max(delta, stats.max_abs_delta);
            stats.max_rel_delta = std::max(delta_rel, stats.max_rel_delta);
#    endif // BENCH_COMPARE_PRINT_STATISTICS

            same = false;
        }
        else {
            //if (print && i < 10)
            //  cout <<std::setprecision(20) << std::scientific << array << "["<<i<<"]:\tsame value ("<<d1[i]<<":"<<d2[i]<<")"<<endl;
            // cout << "same value " << endl;
        }
    }
#    ifdef BENCH_COMPARE_PRINT_STATISTICS
    if (stats.num_failures > 0) {
        stats.mean_abs_delta /= stats.num_failures;
        stats.mean_rel_delta /= stats.num_failures;
        stats.mean_val /= stats.num_failures;
        stats.mean_ref /= stats.num_failures;
    }
#    endif // BENCH_COMPARE_PRINT_STATISTICS


    return same;
}
#endif // BENCHMARKING
