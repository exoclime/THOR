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
#include "debug.h"

#ifdef BENCHMARKING
//#include "grid.h"
//#include "planet.h"
#include "esp.h"
#include "grid.h"



binary_test::binary_test(string output_dir_,
                         string output_base_name_):
    output_dir(output_dir_),
    output_base_name(output_base_name_)
{

}

void binary_test::output_reference(ESP & esp,
                                   const string & iteration,
                                   const string & ref_name) {  
        // open file
        string output_name = output_dir + output_base_name
        + ref_name + "_" + iteration + ".h5";
    
        esp.CopyToHost();
        
        storage s(output_name);
        
        s.append_table(esp.Rho_h,
                       esp.nv*esp.point_num,
                       "Density",
                       "kg/m^3");
        
        s.append_table(esp.pressure_h,
                       esp.nv*esp.point_num,
                       "Pressure",
                       "Pa");
        
        s.append_table(esp.Mh_h,
                       esp.nv*esp.point_num*3,
                       "Horizontal Momentum",
                       "kg m/s");
        
        s.append_table(esp.Wh_h,
                       esp.nvi*esp.point_num,
                       "Vertical Momentum",
                       "kg m/s");
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

            
        esp.CopyToHost();

        storage s(output_name, true);

        //cout << "opening " << output_name << endl;
        
        bool density_comp = compare_to_saved_data(s,
                                                  "Density",
                                                  esp.Rho_h,
                                                  esp.nv*esp.point_num);
        
            
        bool pressure_comp = compare_to_saved_data(s,
                                                  "Pressure",
                                                   esp.pressure_h,
                                                   esp.nv*esp.point_num);

        bool h_momentum_comp = compare_to_saved_data(s,
                                                     "Horizontal Momentum",
                                                     esp.Mh_h,
                                                     esp.nv*esp.point_num*3);
         

        bool v_momentum_comp = compare_to_saved_data(s,
                                                     "Vertical Momentum",
                                                     esp.Wh_h,
                                                     esp.nvi*esp.point_num);

        if (! (density_comp && pressure_comp && h_momentum_comp && v_momentum_comp ))
            cout << iteration << "\tref:\t" << ref_name 
                 << "\trho: " << density_comp
                 << "\tP: " << pressure_comp
                 << "\tM_h: " << h_momentum_comp
                 << "\tM_v: " << v_momentum_comp << endl;
        
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
