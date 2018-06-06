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
// Authors: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//          Urs Schroffenegger, Russel Deitrick
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


// Argument types:
// single on switch: e.g. (--dump, --compare?)
// Boolean: HyDiff, rest ?
// int: e.g. num. iteration, GPU_ID
// enum: e.g. HStest?
// string: e.g. initial conditions file, output directory

// do we need positional arguments?

// could have config file as positional argument
// start from rest as default
// with --initial <initial_file> with *.h5 file, start from there as
// initial condition instead of rest
// if no positional argument, use default values?



#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "parser_helpers.h"

// abstract argument interface class
class arg_interface
{
public:

    virtual int get_nargs() = 0;

    virtual int parse_narg(const string & str) = 0;

    virtual bool is_key(const string & str) = 0;
    
    
};

// interface class 
template<typename T>
class arg
{
public:
    arg( const string & short_form_,
         const string & long_form_,
         const T & value_):
        short_form(short_form_),
        long_form(long_form_),
        value(value_),
        is_set(false)
    {
    };

    void parse_narg(const string & str)
    {

    }
    

    int get_nargs()
    {
        return 0;
    }

    bool is_key(const string & str)
    {
        if (short_form == str)
            return true;
        else if (long_form == str)
            return true;
        else
            return false;
    }
    

private:
    string short_form;
    string long_form;

    T value;

    // TODO: check, do we need this? Or always use default?
    // might be needed to know if we want to override values from
    // config file.
    bool is_set;
    
};


// command arguments handling class
class cmdargs
{
public:
    cmdargs();

    // setup by adding arguments
    template<typename T>
    bool add_arg(const string & short_form,
                 const string & long_form,
                 const T & value );

    // parsing arguments
    bool parse(int argc, char ** argv);

    // find argument from long key
    arg_interface * operator[](string idx);
    
        
private:
    // Storage map as key value pair
    vector<std::unique_ptr<arg_interface>> args;

    enum e_parser_state {
        PARSE_FOR_KEY,
        GOT_KEY
    };

    e_parser_state parser_state = PARSE_FOR_KEY;
    
    int parser_state_nargs = 0;
    vector<std::unique_ptr<arg_interface>>::iterator current_arg_interface;

    bool parser_state_machine(const string & str);

    // search predicates for keys in vectors
    //string search_key;
    /*
    bool long_pred(const std::unique_ptr<arg_interface> &a);
    bool short_pred(const std::unique_ptr<arg_interface> &a);
    */
};

// setup by adding arguments
template<typename T>
bool cmdargs::add_arg(const string & short_form,
                      const string & long_form,
                      const T & value
    )    
{
    

    return true;
}

