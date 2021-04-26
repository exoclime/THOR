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
// directories - helper to check and create directories
//
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
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

#pragma once

#include <string>
#include <vector>

using std::string;


std::vector<string> get_files_in_directory(const string& dir_name);

bool match_output_file_numbering_scheme(const string& file_path, string& basename, int& number);

bool create_output_dir(const string& output_dir);

bool path_exists(const string& path);


class path
{
public:
    path(const string& path);

    bool is_absolute() {
        return is_absolute_path;
    }

    // last file extension separated by a '.'
    string suffix();
    // vector of file extensions separated by a '.'
    std::vector<string> suffixes();
    // parts
    std::vector<string> parts();

    // final part of the path
    string name();
    // final part of path, without suffix
    string stem();

    // parent of last element
    string parent();

    string      to_string();
    const char* c_str();


    path& operator/=(const string& rhs) // compound assignment (does not need to be a member,
    {                                   // but often is, to modify the private members)
        elements.push_back(rhs);

        return *this; // return the result by reference
    }

    // friends defined inside class body are inline and are hidden from non-ADL lookup
    friend path operator/(path          lhs, // passing lhs by value helps optimize chained a+b+c
                          const string& rhs) // otherwise, both parameters may be const references
    {
        lhs /= rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }

private:
    string              element_name;
    std::vector<string> elements;

    bool is_absolute_path = false;
};
