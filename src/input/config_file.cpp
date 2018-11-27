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
// Description: basic config file parser
//
//
// Method:
//
// Known limitations: None.
//
//
// Known issues: None.
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffinegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#include "config_file.h"

#include <fstream>
#include <iomanip>
using namespace std;

#include "parser_helpers.h"


config_file::config_file() {
    append_config_var("config_version", version, -1);
}

bool config_file::parse_file(const string& filename) {
    ifstream file(filename, std::ifstream::in);

    if (file.is_open()) {
        return parse_config(file);
    }
    else {
        file.close();
        return false;
    }
}


bool config_file::parse_config(std::basic_istream<char>& config) {

    // initialise key value pairs
    for (auto&& it = config_vars.begin(); it != config_vars.end(); ++it) {
        it->second->set_default();
        // cout << it->first  << "\t" << it->second->to_str() << endl;
    }


    // read input line by line

    bool parse_success = true;


    int cnt = 1;
    for (std::string line; std::getline(config, line);) {

        // cout << left << "line("<<setw(4)<<right<<cnt<<"): [" << line << "]" << setw(45 - line.length()) << " ";


        string key;
        string value;


        if (is_empty_line(line)) {
            // cout << "empty line" << endl;
        }
        else if (is_comment_line(line)) {
            // cout << "comment"<<endl;
        }
        else if (is_key_value_pair(line, key, value)) {
            // cout << "key value pair matched: (" << match[1] << ":" << match[2] <<")"<< endl;
            auto&& it = config_vars.find(key);
            if (it != config_vars.end()) {
                bool parsed = it->second->parse(value);
                if (!parsed) {
                    cout << "ERROR: parsing of value [" << value << "] failed for key [" << key << "] " << endl;
                    parse_success &= false;
                }
            }
            else {
                cout << "WARNING: config file key [" << key << "] does not exist, skipping" << endl;
                // Do not fail if key does not exist
                //parse_success &= false;
            }
        }
        else {
            cout << "ERROR in config file at line " << cnt << " :\t [" << line << "]" << endl;
            parse_success &= false;
        }
        cnt++;
    }

    return parse_success;
}

bool config_file::append_config_var(const string&                           name,
                                    std::unique_ptr<config_entry_interface> entry) {
    // check if entry exists
    auto&& it = config_vars.find(name);
    if (it != config_vars.end()) {
        cout << "config entry " << name << " already exists" << endl;
        return false;
    }
    else {
        // append entry
        config_vars[name] = std::move(entry);
        return true;
    }
}
