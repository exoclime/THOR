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
// Description: test for storage class
//
//
// Method: writes an array to a file and reloads it
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



#include <iostream>
#include <sstream>

#include "config_file.h"
#include "testing.h"

using namespace std;


    
    

int main()
{
    cout << "config file test" << endl;

    config_file cf;

    string config = "test = 0\n"
        "# comment line\n"
        "other = TRUE\n"
        " # moar comments\n"
        "othertest = icle\n"
        "   \t  \n"
        "\t# undefined, check default value\n";
    
        
    cout << config << endl;

    int test = -1;
    config_entry<int> c(test, 2);
    
//    cf.append_config_var<int>("test", c);
    
    cf.parse_config(std::basic_istringstream<char>(config));

    if (test_val<int>("test var", test, 0) )
        cout << "Test SUCCESS" << endl;
    else
        cout << "Test FAIL" << endl;
            

    exit(0);
    
}
