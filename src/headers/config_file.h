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

#include <iostream>
using namespace std;

using std::string;
using std::map;


class config_entry_interface
{
public:
    virtual void set();
    
};

template <typename T>
class config_entry : public config_entry_interface {
public:
    config_entry(T & target_, T default_val_):
        target(target_),
        default_val(default_val_)
    {
        
    }
    
    void set() 
    {

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
    
        
    template<typename T>
    void append_config_var(const string & name,
                           config_entry_interface entry);

    bool parse_config(istringstream config);
private:
    map<string, config_entry_interface> config_vars;
    
    
};

    
template<typename T>
void config_file::append_config_var(const string & name,
                                    config_entry_interface entry)
{
    std::map<string,config_entry_interface>::iterator it;

    // check if entry exists
    it = config_vars.find(name);
    if (it != config_vars.end())
        // append entry
        config_vars[name] = entry;
    else
        cout << "entry " << name << " already exists"<< endl;
}


