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
// Description: test for grid setup
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

#include "grid.h"

#include <iostream>
using namespace std;


int main ()
{
    bool output = true;
    bool compare = false;
    
    Icogrid Grid(true         , // Spring dynamics option
                 1.15         , // Parameter beta for spring dynamics 
                 4            , // Horizontal resolution level
                 32           , // Number of vertical layers
                 5            ,
                 6371000.0    , // Planet radius
                 36000.0      , // Top of the model's domain
                 false);        // sponge layer  
     

    binary_test & btester = binary_test::get_instance();
    btester.set_output("grid_test", "./results/grid_test/");
    
    if (output)
        btester.output_reference_grid(Grid);
    
    if (compare)
    {
        
        bool result = btester.compare_to_reference_grid(Grid);
        if (result)
            cout << "grid compare SUCCESS" << endl;
        else
            cout << "grid compare FAIL" << endl;            
    }
    
    
    exit(0);
}
