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
// Description: test for command line args
//
//
// Method: reads command line args and does stuff in return. Needs
// external calling
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

#include "cmdargs.h"
#include "testing.h"

using namespace std;


    
    

int main(int argc,  char** argv)
{
    cout << "commandline file test" << endl;
    
    cmdargs argparser;
    argparser.add_arg("t", "test", 0);

    argparser.parse(argc, argv);

    //std::unique_ptr<arg_interface> test_arg = argparser["test"];

    //printf("test %d\n", test_arg->value);
    

    bool success = true;
    /*
    success = success & test_val<int>("int1", int1, 0);
        
    if (success)
    {
        cout << "Test PASS" << endl;
        exit(-1);
    }
    else
    {
        cout << "Test FAIL" << endl;
        exit(-1);
    }    
    */
}
