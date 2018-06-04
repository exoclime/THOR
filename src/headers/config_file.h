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
// Description: basic config file class
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
#pragma once

#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <regex>

using namespace std;

using std::string;
using std::map;


class config_entry_interface
{
public:
    virtual bool parse(string value) = 0;

    virtual void set_default() = 0;
    
};


template< typename T>
bool parse_data(const string & value, T & target)
{

    return false;
}


template<>
inline bool parse_data(const string & value, bool & target )
{
    regex bool_regex("^(true|false)$");

    std::smatch match;

    if (regex_match(value, match, bool_regex))
    {
        target = (match[1] == "true");
        // cout << "parsed [" << value << "] as [" << target << "]" << endl;

        return true;
    }

    return false;
}

template<>
inline bool parse_data(const string & value, int & target )
{
    regex int_regex("^((-|\\+)?[0-9]+)$");

    std::smatch match;

    if (regex_match(value, match, int_regex))
    {
        target = std::stoi(match[1]);
        // cout << "parsed [" << value << "] as [" << target << "]" << endl;

        return true;
    }

    return false;
}
    
template<>
inline bool parse_data(const string & value, double & target )
{

    regex double_regex("^((-|\\+)?[0-9]+(\\.[0-9]+)?((E|e)(-|\\+)?[0-9]+)?)$");



    std::smatch match;

    if (regex_match(value, match, double_regex))
    {
        target = std::stod(match[1]);
        // cout << "parsed [" << value << "] as [" << target << "]" << endl;

        return true;
    }

    return false;
}
    

template<>
inline bool parse_data(const string & value, string & target )
{
    target = value;
    
    return true;
}
    

template <typename T>
class config_entry : public config_entry_interface {
public:
    config_entry(T & target_, T default_val_):
        target(target_),
        default_val(default_val_)
    {
        
    }
    
    bool parse(string value)
    {
        if (parse_data(value, target))
        {
            return true;
        }
        else
        {
            target = default_val;
            
            return false;
        }
    }

    void set_default()
    {
        target = default_val;
    }
    
    
        
private:
    T & target;
    T default_val;
    
};

class config_file
{
public:
    config_file();
    enum var_type { INT, DOUBLE, BOOL, FLOAT, STRING };
    
    bool append_config_var(const string & name,
                           std::unique_ptr<config_entry_interface> entry);


    template<typename T>
    bool append_config_var(const string & name, T & target_, const T & default_val_);
    
    bool parse_config(istringstream config);
private:
    map<string, std::unique_ptr<config_entry_interface>> config_vars;
    
    
};

template<typename T>
bool config_file::append_config_var(const string & name,
                                    T & target_,
                                    const T & default_val_)
{
    append_config_var(name, std::unique_ptr<config_entry<T>>(new config_entry<T>(target_, default_val_)));
}

 
