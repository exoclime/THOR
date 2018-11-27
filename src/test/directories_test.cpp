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
// Description: test for directories class
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
#include <vector>

#include "directories.h"
#include "testing.h"


using namespace std;

template<typename T>
bool print_compare(const string&    msg,
                   const vector<T>& val,
                   const vector<T>& ref) {
    bool test = true;

    if (ref.size() != val.size()) {
        test = false;
        cout << "\033[1;31m"
                "ERROR"
                "\033[0m: "
             << "wrong array sizes for comparison expected: "
             << ref.size() << " got: " << val.size() << endl;
    }
    else {
        for (int i = 0; i < ref.size(); i++) {
            if (ref[i] != val[i]) {
                cout << "\033[1;31m"
                        "ERROR"
                        "\033[0m: "
                     << "wrong value for element[" << i
                     << "] expected: " << ref[i]
                     << "\t got: " << val[i] << endl;
                test = false;
            }
        }
    }

    cout << "[" << ((test) ? "\033[1;32m"
                             "PASS"
                             "\033[0m"
                           : "["
                             "\033[1;31m"
                             "FAIL"
                             "\033[0m")
         << "] " << msg << endl;

    return test;
}


bool print_compare(const string& msg,
                   const string& val,
                   const string& ref) {
    bool test = ref == val;

    if (!test) {
        cout << "\033[1;31m"
                "ERROR"
                "\033[0m: "
             << "wrong value"
             << " expected: " << ref
             << "\t got: " << val << endl;
    }

    cout << "[" << ((test) ? "\033[1;32m"
                             "PASS"
                             "\033[0m"
                           : "["
                             "\033[1;31m"
                             "FAIL"
                             "\033[0m")
         << "] " << msg << endl;

    return test;
}


template<typename T>
bool print_compare(const string& msg,
                   const T&      val,
                   const T&      ref) {
    bool test = ref == val;

    if (!test) {
        cout << "\033[1;31m"
                "ERROR"
                "\033[0m: "
             << "wrong value"
             << " expected: " << ref
             << "\t got: " << val << endl;
    }

    cout << "[" << ((test) ? "\033[1;32m"
                             "PASS"
                             "\033[0m"
                           : "["
                             "\033[1;31m"
                             "FAIL"
                             "\033[0m")
         << "] " << msg << endl;

    return test;
}


int main() {
    bool success = true;
    // ***************************************************************
    cout << "path file test" << endl;

    bool path_test = true;

    {
        // absolute
        path p("/usr/share/docs/");

        success &= print_compare("/usr/share/docs/"
                                 " | absolute path",
                                 p.is_absolute(),
                                 true);

        vector<string> parts = p.parts();
        success &= print_compare("/usr/share/docs/"
                                 " | parts",
                                 parts,
                                 {"usr", "share", "docs"});
        success &= print_compare("/usr/share/docs/"
                                 " | name",
                                 p.name(),
                                 "docs");
        success &= print_compare("/usr/share/docs/"
                                 " | parent",
                                 p.parent(),
                                 "/usr/share");
        success &= print_compare("/usr/share/docs/"
                                 " | to_string",
                                 p.to_string(),
                                 "/usr/share/docs");
        p /= "test";
        success &= print_compare("/usr/share/docs/ + test"
                                 " | parts",
                                 p.parts(),
                                 {"usr", "share", "docs", "test"});
        path p2 = p / "subtest";
        success &= print_compare("/usr/share/docs/test/ + subtest"
                                 " | parts",
                                 p2.parts(),
                                 {"usr", "share", "docs", "test", "subtest"});
        // test integrity of original
        success &= print_compare("/usr/share/docs/test"
                                 " | parts",
                                 p.parts(),
                                 {"usr", "share", "docs", "test"});
    }

    {
        // absolute
        path p("/usr/share/docs");

        success &= print_compare("/usr/share/docs"
                                 " | absolute path",
                                 p.is_absolute(),
                                 true);

        vector<string> parts = p.parts();
        success &= print_compare("/usr/share/docs"
                                 " | parts",
                                 parts,
                                 {"usr", "share", "docs"});
        success &= print_compare("/usr/share/docs"
                                 " | name",
                                 p.name(),
                                 "docs");
        success &= print_compare("/usr/share/docs"
                                 " | parent",
                                 p.parent(),
                                 "/usr/share");
        success &= print_compare("/usr/share/docs"
                                 " | to_string",
                                 p.to_string(),
                                 "/usr/share/docs");
    }

    {
        // relative
        path p("temp/bin.exe");

        success &= print_compare("temp/bin.exe"
                                 " | relative path",
                                 p.is_absolute(),
                                 false);
        vector<string> parts = p.parts();
        success &= print_compare("temp/bin.exe"
                                 " | parts",
                                 parts,
                                 {"temp", "bin.exe"});
        success &= print_compare("temp/bin.exe"
                                 " | name",
                                 p.name(),
                                 "bin.exe");
        success &= print_compare("temp/bin.exe"
                                 " | parent",
                                 p.parent(),
                                 "temp");
        success &= print_compare("temp/bin.exe"
                                 " | to_string",
                                 p.to_string(),
                                 "temp/bin.exe");
    }

    {
        // relative
        path p("./temp/bin.exe");

        success &= print_compare("./temp/bin.exe"
                                 " | relative path with .",
                                 p.is_absolute(),
                                 false);
        vector<string> parts = p.parts();
        success &= print_compare("./temp/bin.exe"
                                 " | parts",
                                 parts,
                                 {".", "temp", "bin.exe"});
        success &= print_compare("./temp/bin.exe"
                                 " | name",
                                 p.name(),
                                 "bin.exe");
        success &= print_compare("./temp/bin.exe"
                                 " | parent",
                                 p.parent(),
                                 "./temp");
        success &= print_compare("./temp/bin.exe"
                                 " | to_string",
                                 p.to_string(),
                                 "./temp/bin.exe");
    }

    {
        path p("temp/subfolder/myfile.tar.gz");

        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | name",
                                 p.name(),
                                 "myfile.tar.gz");
        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | suffix",
                                 p.suffix(),
                                 "gz");
        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | stem",
                                 p.stem(),
                                 "myfile.tar");

        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | parent",
                                 p.parent(),
                                 "temp/subfolder");
        vector<string> parts = p.parts();
        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | parts",
                                 parts,
                                 {"temp", "subfolder", "myfile.tar.gz"});


        vector<string> suffixes = p.suffixes();
        success &= print_compare("temp/subfolder/myfile.tar.gz"
                                 " | suffixes",
                                 suffixes,
                                 {"tar", "gz"});
    }


    {
        path p("results/esp_output_Earth_5.h5");


        success &= print_compare("results/esp_output_Earth_5.h5"
                                 " | suffix",
                                 p.suffix(),
                                 "h5");
        success &= print_compare("results/esp_output_Earth_5.h5"
                                 " | name",
                                 p.name(),
                                 "esp_output_Earth_5.h5");
        success &= print_compare("results/esp_output_Earth_5.h5"
                                 " | stem",
                                 p.stem(),
                                 "esp_output_Earth_5");
        success &= print_compare("results/esp_output_Earth_5.h5"
                                 " | parent",
                                 p.parent(),
                                 "results");
    }

    {
        string basename("");
        int    number = -1;


        bool match = match_output_file_numbering_scheme("esp_output_Earth_5.h5",
                                                        basename,
                                                        number);

        success &= print_compare("match check", match, true);
        success &= print_compare("match name", basename, "Earth");
        success &= print_compare("match number", number, 5);
    }

    {
        string basename("");
        int    number = -1;


        bool match = match_output_file_numbering_scheme("esp_output_super_hot_planet_42.h5",
                                                        basename,
                                                        number);

        success &= print_compare("match check", match, true);
        success &= print_compare("match name", basename, "super_hot_planet");
        success &= print_compare("match number", number, 42);
    }

    {
        string basename("");
        int    number = -1;


        bool match = match_output_file_numbering_scheme("esp_ofvgaput_Earth.h5",
                                                        basename,
                                                        number);

        success &= print_compare("match check", match, false);
    }

    {
        string basename("");
        int    number = -1;


        bool match = match_output_file_numbering_scheme("esp_output_Earth_ret.h5",
                                                        basename,
                                                        number);

        success &= print_compare("match check", match, false);
    }


    // ***************************************************************
    cout << "directories file test" << endl;

    bool directory_success = false;

    directory_success = create_output_dir("test_dir/subtest_dir/");
    if (!directory_success) {
        success = false;
        cout << "directory creation failed" << endl;
    }
    // ***************************************************************


    if (success)
        cout << "Test PASS" << endl;
    else
        cout << "Test FAIL" << endl;


    exit(0);
}
