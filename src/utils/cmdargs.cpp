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
// Description: command line args parsing class
//
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
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////


#include "cmdargs.h"


#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <cstdio>

#include "parser_helpers.h"

using std::cout;
using std::endl;


// Bool argument is just a switch
// Number of arguments following this key is 0
template<> int arg<bool>::get_nargs() {
    return 0;
}

// Flag it as set, it will use its predefined value
template<> void arg<bool>::set_default() {
    has_value = true;
}


cmdargs::cmdargs(const string& app_name_, const string& app_desc_) :
    app_name(app_name_),
    app_desc(app_desc_) {
}

bool cmdargs::parse(int argc, char** argv) {
    bool out = true;

    // initialise parser
    parser_state          = PARSE_FOR_KEY;
    parser_state_nargs    = 0;
    current_arg_interface = args.end();

    for (int i = 1; i < argc; i++) {
        out &= parser_state_machine(argv[i]);
    }

    if (parser_state != PARSE_FOR_KEY) {
        cout << "Error: parsing ended before reading all arguments" << endl;
        out = false;
    }


    return out;
}

// Access operator
arg_interface* cmdargs::operator[](string idx) {
    const string s = idx;

    auto&& it =
        std::find_if(args.begin(),
                     args.end(),
                     [s](const std::unique_ptr<arg_interface>& m) -> bool { return m->is_key(s); });


    if (it != args.end()) {
        return it->get();
    }
    else {
        char error[256];
        sprintf(
            error, "%s:(%d) key [%s] not found in argument list", __FILE__, __LINE__, idx.c_str());

        throw std::runtime_error(error);
    }
}

void cmdargs::print_help() {
    cout << "Usage: " << app_name << " [OPTION]... [CONFIGFILE]" << endl;
    cout << app_desc << endl;

    cout << endl;
    cout << "Positional argument:" << endl;
    positional_arg->print_desc();
    cout << endl;

    cout << "Keyword arguments:" << endl;

    for (auto&& a : args)
        a->print_desc();
}


// Main state machine to parse arguments
bool cmdargs::parser_state_machine(const string& str) {
    switch (parser_state) {
        case PARSE_FOR_KEY: {
            // Check argument key
            current_arg_interface = args.end();
            bool match_key        = false;

            string ma = "";

            if (is_short_cmdargs(str, ma) || is_long_cmdargs(str, ma)) {
                if (ma == "h" || ma == "help") {
                    print_help();
                    exit(0);
                }

                current_arg_interface =
                    std::find_if(args.begin(),
                                 args.end(),
                                 [ma](const std::unique_ptr<arg_interface>& m) -> bool {
                                     return m->is_key(ma);
                                 });
                match_key = true;
            }

            // act for argument
            if (match_key) {
                // is a keyword argument, parse its value

                if (current_arg_interface != args.end()) {
                    // Check if we need follow up arguments
                    parser_state_nargs = (*current_arg_interface)->get_nargs();
                    if (parser_state_nargs == 0) {
                        // No follow up argument, set key to argument default value
                        (*current_arg_interface)->set_default();

                        parser_state = PARSE_FOR_KEY;
                    }
                    else {
                        parser_state = GOT_KEY;
                    }


                    return true;
                }
                else {
                    cout << "command line: error finding key: " << str << endl;

                    // not a key, return error
                    return false;
                }
            }
            else {
                // positional argument, append to list of arguments

                bool parsed = positional_arg->parse_narg(str);
                if (!parsed) {
                    cout << "error parsing value " << str << " for positional argument" << endl;
                }
                return parsed;
            }


            break;
        }
        case GOT_KEY: {
            bool parsed = (*current_arg_interface)->parse_narg(str);

            if (!parsed) {
                cout << "error parsing value " << str << " for key "
                     << (*current_arg_interface)->get_name() << endl;
            }

            string ma = "";
            if (is_short_cmdargs(str, ma) || is_long_cmdargs(str, ma)) {
                cout << "received a command switch when expecting an argument" << endl;
                parsed = false;
            }


            parser_state_nargs--;

            if (parser_state_nargs == 0)
                parser_state = PARSE_FOR_KEY;

            return parsed;


            break;
        }
    }


    return false;
}
