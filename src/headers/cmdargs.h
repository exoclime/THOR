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

    virtual bool parse_narg(const string & str) = 0;

    virtual bool is_key(const string & str) = 0;
    virtual string get_name() = 0;    
};

// interface class 
template<typename T>
class arg : public arg_interface
{
public:
    arg( const string & short_form_,
         const string & long_form_,
         const T & value_):
        short_form(short_form_),
        long_form(long_form_),
        value(value_),
        has_value(false)
    {
    };

    bool parse_narg(const string & str)
    {
        if (parse_data(str, value))
        {
            has_value = true;
            return true;
        }
        else
        {
            has_value = false;
            return false;
        }
        
    }
    

    int get_nargs()
    {
        return 1;
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

    T get()
    {
        return value;
    }

    bool is_set()
    {
        return has_value;
        
    }

    string get_name() 
    {
        return long_form;
    }
    
    
    

private:
    string short_form;
    string long_form;

    T value;

    // TODO: check, do we need this? Or always use default?
    // might be needed to know if we want to override values from
    // config file.
    bool has_value;
    
    
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

    template<typename T>
    bool get_arg(const string & long_form,
                 T & value );

    // positional argument
    template<typename T>
    bool add_positional_arg(const T & value );

    template<typename T>
    bool get_positional_arg( T & value );

    // parsing arguments
    bool parse(int argc, char ** argv);
    
        
private:
     // find argument from long key
    arg_interface * operator[](string idx);
    
    // Storage map as key value pair
    vector<std::unique_ptr<arg_interface>> args;

    // Storage for positional argument
    std::unique_ptr<arg_interface> positional_arg = nullptr;
    
    enum e_parser_state {
        PARSE_FOR_KEY,
        GOT_KEY
    };

    e_parser_state parser_state = PARSE_FOR_KEY;
    
    int parser_state_nargs = 0;
    vector<std::unique_ptr<arg_interface>>::iterator current_arg_interface;

    bool parser_state_machine(const string & str);

    void add_arg(std::unique_ptr<arg_interface> arg)
    {
        args.push_back(std::move(arg));
    }
    
};

// setup by adding argumentspositional argument
template<typename T>
bool cmdargs::add_positional_arg(
    const T & value
    )    
{
    positional_arg = std::unique_ptr<arg<T>>(new arg<T>("", "", value));
    
    return true;
}

template<typename T>
bool cmdargs::get_positional_arg( T & value    )    
{

    if (positional_arg != nullptr)
    {
        arg<T> * argument  = (arg<T>* )(positional_arg.get());

        value = argument->get();

        return argument->is_set();
    }
    else
    {
        throw std::runtime_error( "no positional argument defined" );
    }
}

// setup by adding arguments
template<typename T>
bool cmdargs::add_arg(const string & short_form,
                      const string & long_form,
                      const T & value
    )    
{
    // check that the key wasn't already added

    // check short form

    auto && it = std::find_if(args.begin(),
                              args.end(),
                              [short_form](const std::unique_ptr<arg_interface>  & m)
                              -> bool { return m->is_key(short_form); });
    
                              
    if (it != args.end())
    {
        char error[256];
        sprintf(error,
                "%s:(%d) key %s already exists in argument parser",
                __FILE__,
                __LINE__,
                short_form);
        
        throw std::range_error(error );        
    }
    auto && it2 = std::find_if(args.begin(),
                               args.end(),
                               [long_form](const std::unique_ptr<arg_interface>  & m)
                               -> bool { return m->is_key(long_form); });
    if (it2 != args.end())
    {
        char error[256];
        sprintf(error,
                "%s:(%d) key %s already exists in argument parser",
                __FILE__,
                __LINE__,
                long_form);
        
        throw std::range_error(error );        
    }

    // add the new key
    add_arg(std::unique_ptr<arg<T>>(new arg<T>(short_form, long_form, value)));
    
    
    return true;
}

// setup by adding arguments
template<typename T>
bool cmdargs::get_arg(const string & long_form,
                      T & value    )    
{
    arg<T> * argument  = (arg<T> *)(*this)[long_form];

    value = argument->get();
    
    

    return argument->is_set();
}

