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

#include <iomanip>
#include <sstream>

binary_test::binary_test(string output_dir_,
                         string output_base_name_):
    output_dir(output_dir_),
    output_base_name(output_base_name_)
{

}

void binary_test::init_output_data(ESP & esp)
{
    if (output_defined)
        return;

    // define all the data we want to output
    output_definitions =
        {
            { esp.Rho_d,         esp.nv*esp.point_num,   "Density", "rho", true },
            { esp.pressure_d,    esp.nv*esp.point_num,   "Pressure", "P", true },
            { esp.Mh_d,          esp.nv*esp.point_num*3, "Horizontal Momentum", "Mh", true},
            { esp.Wh_d,          esp.nvi*esp.point_num,  "Vertical Momentum", "Wh", true},
            { esp.temperature_d, esp.nv*esp.point_num,   "Temperature", "T", true},
            { esp.W_d,           esp.nv*esp.point_num,   "W vert momentum", "W", true },
            // RK variables
            { esp.pressures_d,   esp.nv*esp.point_num,   "RK pressures", "ps", true},
            { esp.Rhos_d,        esp.nv*esp.point_num,   "RK Rhos", "rhos", true},
            { esp.Mhs_d,         esp.nv*esp.point_num*3, "RK Mhs", "Mhs", true},
            { esp.Ws_d,          esp.nv*esp.point_num,   "RK Ws", "Ws", true},
            { esp.Whs_d,         esp.nvi*esp.point_num,  "RK Whs", "Whs", true},
            { esp.pressurek_d,   esp.nv*esp.point_num,   "RK pressurek", "pk", true},
            { esp.Rhok_d,        esp.nv*esp.point_num,   "RK Rhok", "Rhok", true},
            { esp.Mhk_d,         esp.nv*esp.point_num*3, "RK Mhk", "Mhk", true},
            { esp.Wk_d,          esp.nv*esp.point_num,   "RK Wk", "Wk", true},
            { esp.Whk_d,         esp.nvi*esp.point_num,  "RK Whk", "Whk", true},
        };

    // allocate temporary buffer
    int max_size = 0;
    
    for (auto & def : output_definitions)
    {
        if (def.size > max_size)
            max_size = def.size;
    }
    
    mem_buf = std::unique_ptr<double[]>(new double[max_size], std::default_delete<double[]>());
        
    output_defined = true;
}


void binary_test::output_reference(ESP & esp,
                                   const string & iteration,
                                   const string & ref_name) {  
        // open file
        string output_name = output_dir + output_base_name
        + ref_name + "_" + iteration + ".h5";

        init_output_data(esp);
        
        storage s(output_name);

        for (auto & def : output_definitions)
        {
            // copy data to host if needed
            // and write it to the output file
            if (def.device_ptr)
            {
                esp.getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
                s.append_table(mem_buf.get(),
                               def.size,
                               def.name,
                               "-");
            }
            else
                s.append_table(def.data,
                               def.size,
                               def.name,
                               "-");
        }

}

void binary_test::output_reference_grid(Icogrid & grid) {
    // open file
    string output_name = output_dir + output_base_name
        + "_grid.h5";    
        
    storage s(output_name);

    s.append_table(grid.point_xyz,
                   grid.point_num*3,
                   "xyz",
                   "m");
    
    s.append_table(grid.point_xyzq,
                   grid.point_num*3*6,
                   "xyzq",
                   "m");
    
    s.append_table(grid.pent_ind,
                   12,
                   "pent_ind",
                   "-");
    
    s.append_table(grid.point_local,
                   6*grid.point_num,
                   "point_local",
                   "-");
    
    s.append_table(grid.halo,
                   grid.nh,
                   "halo",
                   "-");
    
    s.append_table(grid.maps,
                   (grid.nl_region+2)*(grid.nl_region+2)*grid.nr,
                   "maps",
                   "-");
    
    s.append_table(grid.func_r,
                   3*grid.point_num,
                   "func_r",
                   "-");
    
    s.append_table(grid.areas,
                   6*3*grid.point_num,
                   "areas",
                   "-");
    
    s.append_table(grid.areasTr,
                   6*grid.point_num,
                   "areasTr",
                   "-");
    
    s.append_table(grid.areasT,
                   grid.point_num,
                   "areasT",
                   "-");
    
    s.append_table(grid.nvec,
                   6*3*grid.point_num,
                   "nvec",
                   "-");
    s.append_table(grid.nvecoa,
                   6*3*grid.point_num,
                   "nvecoa",
                   "-");
    s.append_table(grid.nvecti,
                   6*3*grid.point_num,
                   "nvecti",
                   "-");
    s.append_table(grid.nvecte,
                   6*3*grid.point_num,
                   "nvecte",
                   "-");
    
    s.append_table(grid.Altitude,
                   grid.nv,
                   "Altitude",
                   "m");
    s.append_table(grid.Altitudeh,
                       grid.nvi,
                   "Altitudeh",
                   "m");
    
    s.append_table(grid.lonlat,
                   2*grid.point_num,
                   "lonlat",
                       "°");
    s.append_table(grid.div,
                   7*3*grid.point_num,
                   "div",
                       "-");
    s.append_table(grid.grad,
                   7*3*grid.point_num,
                   "grad",
                   "-");    
}

binary_test & binary_test::get_instance() {
    static binary_test bt(
                          BENCHMARK_DUMP_REF_PATH,
                          BENCHMARK_DUMP_BASENAME);

    
    return bt;
}


bool binary_test::compare_to_reference(ESP & esp,
                                       const string & iteration,
                                       const string & ref_name) {    
  
        string output_name = output_dir + output_base_name
        + ref_name + "_" + iteration + ".h5";

        storage s(output_name, true);

        

        init_output_data(esp);
        
        bool out = true;
        std::ostringstream oss;
        oss << std::left << std::setw(50) << output_name;
        oss  << std::setw(6) << iteration << " ref: " << std::setw(30) << ref_name;
        
        for (auto & def : output_definitions)
        {
            bool comp = false;
            
            // copy data to host if needed
            // and write it to the output file
            if (def.device_ptr)
            {
                esp.getDeviceData(def.data, mem_buf.get(), def.size * sizeof(double));
                comp = compare_to_saved_data(s, def.name, mem_buf.get(), def.size);
            }
            else
            {
                comp = compare_to_saved_data(s, def.name, def.data, def.size);
            }

            oss << " " << def.short_name <<": " << comp;
            out &= comp;
            
            
        }
        
                                              
                                              
                                              
            

//        if (!out )
            cout << oss.str() << endl;
        
        return out;
}



bool binary_test::compare_to_reference_grid(Icogrid & grid) {
    string output_name = output_dir + output_base_name
        + "_grid.h5";
    
    storage s(output_name, true);

    //cout << "opening " << output_name << endl;
    
    bool out = true;
    out &= compare_to_saved_data(s,
                                 "xyz",
                                 grid.point_xyz,
                                 grid.point_num*3);
    
    out &= compare_to_saved_data(s,
                                 "xyzq",
                                 grid.point_xyzq,
                                 grid.point_num*3*6);
    

    
    out &= compare_to_saved_data(s,
                                 "pent_ind",
                                 grid.pent_ind,
                                 12);
    
    out &= compare_to_saved_data(s,
                                 "point_local",
                                 grid.point_local,
                                 6*grid.point_num);

    out &= compare_to_saved_data(s,
                                 "halo",
                                 grid.halo,
                                 grid.nh);
    
    out &= compare_to_saved_data(s,
                                 "maps",
                                 grid.maps,
                                 (grid.nl_region+2)*(grid.nl_region+2)*grid.nr);
    
    out &= compare_to_saved_data(s,
                                 "func_r",
                                 grid.func_r,
                                 3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "areas",
                                 grid.areas,
                                 6*3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "areasTr",
                                 grid.areasTr,
                                 6*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "areasT",
                                 grid.areasT,
                                 grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "nvec",
                                 grid.nvec,
                                 6*3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "nvecoa",
                                 grid.nvecoa,
                                 6*3*grid.point_num );
    
    out &= compare_to_saved_data(s,
                                 "nvecti",
                                 grid.nvecti,
                                 6*3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "nvecte",
                                 grid.nvecte,
                                 6*3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "Altitude",
                                 grid.Altitude,
                                 grid.nv);
    out &= compare_to_saved_data(s,
                                 "Altitudeh",
                                 grid.Altitudeh,
                                 grid.nvi);
    
    out &= compare_to_saved_data(s,
                                 "lonlat",
                                 grid.lonlat,
                                 2*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "div",
                                 grid.div,
                                 7*3*grid.point_num);
    
    out &= compare_to_saved_data(s,
                                 "grad",
                                 grid.grad,
                                 7*3*grid.point_num);
    
    return out;    
}

#endif // BENCHMARKING
