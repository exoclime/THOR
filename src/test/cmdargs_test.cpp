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


int main(int argc, char** argv) {
    cout << "commandline file test" << endl;

    cmdargs argparser("cmdargs_test", "tester for argument parser");

    argparser.add_positional_arg(string("positional"), "positional argument");

    argparser.add_arg("i", "int", -1, "int argument");
    argparser.add_arg("b", "bool", true, "boolean argument");
    argparser.add_arg("n", "negbool", false, "negative boolean argument");
    argparser.add_arg("d", "double", 1e-5, "double argument");
    argparser.add_arg("s", "string", string("test"), "string argument");
    argparser.parse(argc, argv);

    // positional args
    string positional_arg_string     = "fart";
    bool   positional_arg_string_set = argparser.get_positional_arg(positional_arg_string);
    printf("# positional argument - value: %s\tset: %d\n", positional_arg_string.c_str(), positional_arg_string_set);


    // keyqord args
    int  test_arg_int     = -100;
    bool test_arg_int_set = argparser.get_arg("int", test_arg_int);
    printf("int - value: %d\tset: %d\n", test_arg_int, test_arg_int_set);

    // In this example case, boolean is false by default and key puts it to true
    bool test_arg_bool     = false;
    bool test_arg_bool_set = argparser.get_arg("bool", test_arg_bool);
    printf("bool - value: %d\tset: %d\n", test_arg_bool, test_arg_bool_set);

    // In this example case, boolean is true by default and key puts it to false
    bool test_arg_negbool     = true;
    bool test_arg_negbool_set = argparser.get_arg("negbool", test_arg_negbool);
    printf("negbool - value: %d\tset: %d\n", test_arg_negbool, test_arg_negbool_set);

    double test_arg_double     = -1e5;
    bool   test_arg_double_set = argparser.get_arg("double", test_arg_double);
    printf("double - value: %e\tset: %d\n", test_arg_double, test_arg_double_set);

    string test_arg_string     = "fart";
    bool   test_arg_string_set = argparser.get_arg("string", test_arg_string);
    printf("string - value: %s\tset: %d\n", test_arg_string.c_str(), test_arg_string_set);
}
