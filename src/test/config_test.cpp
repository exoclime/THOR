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


#include <iomanip>
#include <iostream>
#include <sstream>

#include "config_file.h"
#include "testing.h"

using namespace std;


int main() {
    cout << "config file test" << endl;

    config_file cf;

    string config = "int1 = 0\n"
                    " int2=  1234567890 \n"
                    "int3\t=-1234567890\t\n"
                    "# comment line\n"
                    "other = false\n"
                    " # moar comments\n"
                    "othertest = icle\n"
                    "onemore = true # line with comment\n"
                    "   \t  \n"
                    "\t# undefined, check default value\n"
                    "double1 = 1.0\n"
                    "double2 = -3.14e19\n"
                    "double3 = 1242E+32\n";

    int           i = 1;
    istringstream iss(config);


    for (string line; getline(iss, line);) {

        cout << setw(4) << i << ":\t" << line << endl;
        i++;
    }


    cout << "appending configurations" << endl;
    int int1 = -1;
    cf.append_config_var("int1", int1, 2);
    int int2 = -1;
    cf.append_config_var("int2", int2, 2);
    int int3 = -1;
    cf.append_config_var("int3", int3, 2);

    double double1 = 0.0;
    cf.append_config_var("double1", double1, -1.0);
    double double2 = 0.0;
    cf.append_config_var("double2", double2, -1.0);
    double double3 = 0.0;
    cf.append_config_var("double3", double3, -1.0);

    bool onemore = false;
    cf.append_config_var("onemore", onemore, false);
    bool other = true;
    cf.append_config_var("other", other, true);
    string othertest("fart");
    cf.append_config_var("othertest", othertest, string("tickle"));

    // default values
    bool def_bool = false;
    cf.append_config_var("def_bool", def_bool, true);
    int def_int = -2;
    cf.append_config_var("def_int", def_int, 42);
    double def_double = 3.14;
    cf.append_config_var("def_double", def_double, 2.18);
    string def_string("wrong");
    cf.append_config_var("def_string", def_string, string("default"));


    cout << "start parsing" << endl << endl;

    std::basic_istringstream<char> istr(config);


    bool success = true;
    success &= cf.parse_config(istr);
    success &= test_val<int>("int1", int1, 0);
    success &= test_val<int>("int2", int2, 1234567890);
    success &= test_val<int>("int3", int3, -1234567890);
    success &= test_val<double>("double1", double1, 1.0);
    success &= test_val<double>("double2", double2, -3.14e19);
    success &= test_val<double>("double3", double3, 1242e32);
    success &= test_val<bool>("onemore", onemore, true);
    success &= test_val<bool>("other", other, false);
    success &= test_val<string>("othertest", othertest, string("icle"));
    success &= test_val<bool>("def_bool", def_bool, true);
    success &= test_val<int>("def_int", def_int, 42);
    success &= test_val<double>("def_double", def_double, 2.18);
    success &= test_val<string>("def_string", def_string, string("default"));

    if (success)
        cout << "Test PASS" << endl;
    else
        cout << "Test FAIL" << endl;


    exit(0);
}
