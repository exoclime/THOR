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

#include "../headers/binary_test.h"

//#include "grid.h"
//#include "planet.h"
#include "../headers/esp.h"



binary_test::binary_test(ESP & esp_,
                         string output_dir_,
                         string output_base_name_):
    esp(esp_),
    output_dir(output_dir_),
    output_base_name(output_base_name_)
{

}

void binary_test::output_reference(const int & iteration,
                                   const string & ref_name,
                                   const binary_test_mode & mode) {
// open file
    string output_name = output_dir + output_base_name
        + ref_name + "_" + std::to_string(iteration);
    
    if (mode == binary_test_mode::data)
    {
        

    }
    
}

binary_test & binary_test::get_instance(ESP & esp_) {
    static binary_test bt(esp_, "../ref/", "bindata_");
    
    return bt;
}


bool binary_test::compare_to_reference(const int & iteration,
                                       const string & ref_name,
                                       const binary_test_mode & mode) {
    return true;    
}

