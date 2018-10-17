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
//
//
// Known limitations: None.
//
// Known issues: None.
//   
//
// If you use this code please cite the following reference: 
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016  
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "binary_test.h"

#ifdef BENCHMARKING
//#include "grid.h"
//#include "planet.h"
#include "esp.h"
#include "grid.h"
#include "debug_helpers.h"
#include "directories.h"
#include <iomanip>
#include <sstream>

using namespace std;


// define all the variables we can check in the debugging functions
map<string, output_def> build_definitions( ESP & esp, Icogrid & grid)
{
  
    map<string, output_def> out =
        {
            // {"map name, {variable pointer, table size, name, short name, on device}}
	  {"Rho_d",         { esp.Rho_d,         esp.nv*esp.point_num,   "Density", "rho", true}},
	  
	  {"pressure_d",    { esp.pressure_d,    esp.nv*esp.point_num,   "Pressure", "P", true}},
	  {"Mh_d",          { esp.Mh_d,          esp.nv*esp.point_num*3, "Horizontal Momentum", "Mh", true}},
	  {"Wh_d",          { esp.Wh_d,          esp.nvi*esp.point_num,  "Vertical Momentum", "Wh", true}},
	  {"temperature_d", { esp.temperature_d, esp.nv*esp.point_num,   "Temperature", "T", true}},
	  
	  {"W_d",           { esp.W_d,           esp.nv*esp.point_num,   "W vert momentum", "W", true }},
          
	  {"h_d",           { esp.h_d,           esp.nv*esp.point_num,   "Entalphy h", "h", true }},
          {"hh_d",          { esp.hh_d,          esp.nvi*esp.point_num, "Entalphy hh", "hh", true }},
	  // RK variables
	  {"pressures_d",   { esp.pressures_d,   esp.nv*esp.point_num,   "RK pressures", "ps", true}},
	  {"Rhos_d",        { esp.Rhos_d,        esp.nv*esp.point_num,   "RK Rhos", "rhos", true}},
	  {"Mhs_d",         { esp.Mhs_d,         esp.nv*esp.point_num*3, "RK Mhs", "Mhs", true}},
	  {"Ws_d",          { esp.Ws_d,          esp.nv*esp.point_num,   "RK Ws", "Ws", true}},
	  {"Whs_d",         { esp.Whs_d,         esp.nvi*esp.point_num,  "RK Whs", "Whs", true}},
	  {"pressurek_d",   { esp.pressurek_d,   esp.nv*esp.point_num,   "RK pressurek", "pk", true}},
	  {"Rhok_d",        { esp.Rhok_d,        esp.nv*esp.point_num,   "RK Rhok", "Rhok", true}},
	  
	  {"Mhk_d",         { esp.Mhk_d,         esp.nv*esp.point_num*3, "RK Mhk", "Mhk", true}},
	  {"Wk_d",          { esp.Wk_d,          esp.nv*esp.point_num,   "RK Wk", "Wk", true}},
	  {"Whk_d",         { esp.Whk_d,         esp.nvi*esp.point_num,  "RK Whk", "Whk", true}},
	  // local variables
	  
	  {"Adv_d",         { esp.Adv_d,         esp.nv*esp.point_num*3,   "Advection", "Adv", true}},
	  {"v_d",           { esp.v_d,           esp.nv*esp.point_num*3,   "Velocity", "v", true}},
	  
	  // grid
	  {"point_xyz",     { grid.point_xyz, grid.point_num*3, "xyz", "x", false}},
	  
	  {"point_xyzq",    { grid.point_xyzq, grid.point_num*3*6, "xyzq", "xq", false}},
	  
	  /*
              disable int tables for now
              {"pent_ind",      { grid.pent_ind, 12, "pent_ind", "pi", false}},
            
              {"point_local",   { grid.point_local, 6*grid.point_num, "point_local", "pl", false}},
              {"halo",          { grid.halo, grid.nh, "halo", "halo", false}},
              {"maps",          { grid.maps, (grid.nl_region+2)*(grid.nl_region+2)*grid.nr, "maps", "m", false}},
            */
	  
	  {"func_r",        { grid.func_r, 3*grid.point_num, "func_r", "f", false}},
	  {"areas",         { grid.areas, 6*3*grid.point_num, "areas", "a", false}},
	  {"areasTr",       { grid.areasTr, 6*grid.point_num, "areasTr", "aTr", false}},
	  
	  {"areasT",        { grid.areasT, grid.point_num, "areasT", "aT", false}},
	  {"nvec",          { grid.nvec, 6*3*grid.point_num, "nvec", "nc", false}},
	  {"nvecoa",        { grid.nvecoa, 6*3*grid.point_num, "nvecoa", "na", false}},
	  {"nvecti",        { grid.nvecti, 6*3*grid.point_num, "nvecti", "nti", false}},
	  {"nvecte",        { grid.nvecte, 6*3*grid.point_num, "nvecte", "nte", false}},
	  
	  {"Altitude",      { grid.Altitude, grid.nv, "Altitude", "Alt", false}},
	  {"Altitudeh",     { grid.Altitudeh, grid.nvi, "Altitudeh", "Alth", false}},
	  {"lonlat",        { grid.lonlat, 2*grid.point_num, "lonlat", "ll", false}},
	  
	  {"div",           { grid.div, 7*3*grid.point_num, "div", "d", false}},
	  {"grad",          { grid.grad, 7*3*grid.point_num, "grad", "g", false }},
        };

    return out;        
}



binary_test::binary_test(string output_dir_,
                         string output_base_name_):
    output_dir(output_dir_),
    output_base_name(output_base_name_),
    nan_check_d(nullptr)
{
    create_output_dir(output_dir);
    
}
binary_test::~binary_test()
{
    if (nan_check_d != nullptr)
        deinit_device_mem_check(nan_check_d);
}
// generic data test function
void binary_test::check_data(const string & iteration,
                    const string & ref_name,
                    const vector<string> & input_vars,
                    const vector<string> & output_vars)
{    
#ifdef BENCH_CHECK_LAST_CUDA_ERROR
    check_last_cuda_error(ref_name);
#endif // BENCH_CHECK_LAST_CUDA_ERROR
    
    // load definitions for variables
    vector<output_def> data_output;
    for (auto & name: output_vars)
    {
        auto && it = output_definitions.find(name);
        if (it != output_definitions.end())
        {
            data_output.push_back(it->second);
        }
    }

    
// Specific debugging functions
#ifdef BENCH_POINT_WRITE
    output_reference(iteration, ref_name, data_output);
#endif // BENCH_POINT_WRITE

#ifdef BENCH_POINT_COMPARE
    compare_to_reference(iteration, ref_name, data_output);

#endif // BENCH_POINT_COMPARE

#ifdef BENCH_NAN_CHECK
    if ( nan_check_d == nullptr)
    {
        nan_check_d = init_device_mem_check(nan_check_d);
    }
    
    check_nan(iteration, ref_name, data_output);
#endif // BENCH_NAN_CHECK

}
// Check for NaNs
// Returns true if data is ok (contains no NaNs!)
bool binary_test::check_nan(const string & iteration,
                            const string & ref_name,
                            const vector<output_def> & data_output)
{
    bool out = true;
    std::ostringstream oss;
    oss  << std::setw(6) << iteration << " NaNChk: " << std::setw(30) << ref_name;
    
    for (auto & def : data_output)
    {
        bool nan = false;
        // copy data to host if needed
        // and write it to the output file
        nan = check_array_for_nan(def.data, def.size, def.device_ptr, nan_check_d);
        

        oss << " " << def.short_name <<": " << nan;
        out &= !nan; // need to use binary neg, nan is true if data has NaNs, out should be true if
                     // data is OK
        
    }

#ifndef BENCH_PRINT_DEBUG   // if unset, print only failures
        if (!out )
#endif // BENCH_PRINT_DEBUG
            cout << oss.str() << endl;
        
        return out;
}


void binary_test::output_reference(const string & iteration,
                                   const string & ref_name,
                                   const vector<output_def> & data_output) {  
        // open file
        string output_name = output_dir + output_base_name
        + ref_name + "_" + iteration + ".h5";
        
        storage s(output_name);

        for (auto & def : data_output)
        {
            // copy data to host if needed
            // and write it to the output file
	    if (def.device_ptr)
            {
	      getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
              check_last_cuda_error(string("output ") +ref_name + string(" ") +def.name);
              
              
	       s.append_table(mem_buf.get(),
                              def.size,
                              def.short_name,
                              "-",
                              def.short_name);               
            }
            else
	      s.append_table(def.data,
                             def.size,
                             def.short_name,
                             "-",
                             def.short_name);
        }



}

binary_test & binary_test::get_instance() {
    static binary_test bt(
                          BENCHMARK_DUMP_REF_PATH,
                          BENCHMARK_DUMP_BASENAME);

    
    return bt;
}


bool binary_test::compare_to_reference(const string & iteration,
                                       const string & ref_name,
                                       const vector<output_def> & data_output) {    
  
        string output_name = output_dir + output_base_name
        + ref_name + "_" + iteration + ".h5";

        if (!path_exists(output_name))
        {
            cout << "No compare file: " << output_name << endl;
            return true;

        }
        
        storage s(output_name, true);
        
        bool out = true;
        std::ostringstream oss;
        oss << std::left << std::setw(50) << output_name;
        oss  << std::setw(8) << iteration << " ref: " << std::setw(30) << ref_name;
        
        for (auto & def : data_output)
        {
            bool comp = false;
            
            // copy data to host if needed
            // and write it to the output file
            if (!s.has_table(def.short_name))
            {
                // no table in input
                oss << " " << def.short_name <<": " << "N/A";
                
                continue;
            }
            
            if (def.device_ptr)
            {
	      getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
                comp = compare_to_saved_data(s, def.short_name, mem_buf.get(), def.size);
            }
            else
            {
	      comp = compare_to_saved_data(s, def.short_name, def.data, def.size);
            }

            oss << " " << def.short_name <<": " << comp;
            out &= comp;
            
            
        }
#ifndef BENCH_PRINT_DEBUG   // if unset, print only failures
        if (!out )
#endif
            cout << oss.str() << endl;
        
        return out;
}

void binary_test::set_definitions(const map<string, output_def> & defs)
{
  output_definitions = defs;
  int memsize = 0;
  for (auto & d : defs)
  {
    if (d.second.size > memsize)
      memsize = d.second.size;
  }

  mem_buf = std::unique_ptr<double[]>(new double[memsize], std::default_delete<double[]>());
}


#endif // BENCHMARKING
