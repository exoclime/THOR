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
// Description: basic config file parsing class
//
//
// Method: parses a config file and reads it's value to defined variables
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
#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <string>


#include "parser_helpers.h"


using std::map;
using std::string;

// Base interface class for config storage class
class config_entry_interface
{
public:
    virtual bool parse(string value) = 0;

    virtual void set_default() = 0;

    virtual string to_str() = 0;
};


// config entry class, storing data, default value  and calling the
// parsing function
template<typename T> class config_entry : public config_entry_interface
{
public:
    config_entry(T& target_, T default_val_) : target(target_), default_val(default_val_) {
    }

    bool parse(string value) {
        if (parse_data(value, target)) {
            return true;
        }
        else {
            target = default_val;

            return false;
        }
    }

    void set_default() {
        target = default_val;
    }

    string to_str() {
        return to_strg(target);
    }


private:
    T& target;
    T  default_val;
};

// Config file parsing class ,storing the entries and parsing files
class config_file
{
public:
    config_file();

    // config vars definition functions
    template<typename T>
    bool append_config_var(const string& name, T& target_, const T& default_val_);

    // parsing functions
    bool parse_config(std::basic_istream<char>& config);
    bool parse_file(const string& filename);

private:
    // internal parsing function
    bool append_config_var(const string& name, std::unique_ptr<config_entry_interface> entry);

    // Storage map as key value pair
    map<string, std::unique_ptr<config_entry_interface>> config_vars;

    int version = -1;
};

template<typename T>
bool config_file::append_config_var(const string& name, T& target_, const T& default_val_) {
    return append_config_var(
        name, std::unique_ptr<config_entry<T>>(new config_entry<T>(target_, default_val_)));
}
