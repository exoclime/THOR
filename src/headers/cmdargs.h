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
#pragma once


#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "log_writer.h"
#include "parser_helpers.h"

using std::string;

// abstract argument interface class
class arg_interface
{
public:
    virtual int get_nargs() = 0;

    virtual bool parse_narg(const string& str) = 0;

    virtual bool   is_key(const string& str) = 0;
    virtual string get_name()                = 0;
    virtual void   print_desc()              = 0;

    virtual void set_default() = 0;
};

// base template interface class
template<typename T> class arg : public arg_interface
{
public:
    arg(const string& short_form_,
        const string& long_form_,
        const T& value_, // To give the type of the argument, value doesn't get used, other than
                         // for bool that gets set on command line without argument
        const string& help_string_) :
        short_form(short_form_),
        long_form(long_form_),
        help_string(help_string_),
        value(value_),
        has_value(false){};

    // parse argument string to internal value
    // calls template function to convert to correct value
    bool parse_narg(const string& str) {
        if (parse_data(str, value)) {
            has_value = true;
            return true;
        }
        else {
            has_value = false;
            return false;
        }
    };

    // For arguments that have no value, set as having a value so that it uses its default
    void set_default() {
        has_value = true;
    };


    // Number of arguments following this key
    int get_nargs() {
        return 1;
    }

    // check if this argument matches the key
    bool is_key(const string& str) {
        if (short_form == str)
            return true;
        else if (long_form == str)
            return true;
        else
            return false;
    }

    // return the value
    T get() {
        return value;
    }

    // returns if this argument has been set through command line
    bool is_set() {
        return has_value;
    }

    // name of argument
    string get_name() {
        return long_form;
    }

    // description for help string
    void print_desc() {
        if (short_form == "")
            printf("%28s%s\n", " ", help_string.c_str());
        else
            printf(
                " -%s / --%-20s %s\n", short_form.c_str(), long_form.c_str(), help_string.c_str());
    }


private:
    string short_form;
    string long_form;
    string help_string;

    T value;

    // to test if value has been set
    bool has_value;
};


// command arguments handling class
class cmdargs
{
public:
    // Constructor, takes app name and description string as argument
    // for help string
    cmdargs(const string& app_name, const string& app_desc);

    // setup by adding arguments
    template<typename T>
    bool
    add_arg(const string& short_form, const string& long_form, const T& value, const string& help);

    template<typename T> bool get_arg(const string& long_form, T& value);

    // positional argument
    template<typename T> bool add_positional_arg(const T& value, const string& help);

    template<typename T> bool get_positional_arg(T& value);

    // parse argument array
    bool parse(int argc, char** argv);

    // print out help string
    void print_help();

private:
    // find argument from long key
    arg_interface* operator[](string idx);

    // helper to add argument
    void add_arg(std::unique_ptr<arg_interface> arg) {
        args.push_back(std::move(arg));
    }

    // Storage map as key value pair
    std::vector<std::unique_ptr<arg_interface>> args;

    // Storage for positional argument
    std::unique_ptr<arg_interface> positional_arg = nullptr;

    // Parsing state machine
    enum e_parser_state { PARSE_FOR_KEY, GOT_KEY };

    e_parser_state parser_state = PARSE_FOR_KEY;

    int                                                   parser_state_nargs = 0;
    std::vector<std::unique_ptr<arg_interface>>::iterator current_arg_interface;

    bool parser_state_machine(const string& str);


    string app_name;
    string app_desc;
};

// setup by adding positional argument
template<typename T> bool cmdargs::add_positional_arg(const T& value, const string& help) {
    positional_arg = std::unique_ptr<arg<T>>(new arg<T>("", "", value, help));

    return true;
}

template<typename T> bool cmdargs::get_positional_arg(T& value) {

    if (positional_arg != nullptr) {
        arg<T>* argument = (arg<T>*)(positional_arg.get());

        value = argument->get();

        return argument->is_set();
    }
    else {
        throw std::runtime_error("no positional argument defined");
    }
}

// setup by adding arguments
template<typename T>
bool cmdargs::add_arg(const string& short_form,
                      const string& long_form,
                      const T&      value,
                      const string& help) {
    // check that the key wasn't already added

    // check short form

    auto&& it = std::find_if(
        args.begin(), args.end(), [short_form](const std::unique_ptr<arg_interface>& m) -> bool {
            return m->is_key(short_form);
        });


    if (it != args.end()) {
        char error[256];
        sprintf(error,
                "%s:(%d) key %s already exists in argument parser",
                __FILE__,
                __LINE__,
                short_form.c_str());

        throw std::range_error(error);
    }
    auto&& it2 = std::find_if(
        args.begin(), args.end(), [long_form](const std::unique_ptr<arg_interface>& m) -> bool {
            return m->is_key(long_form);
        });
    if (it2 != args.end()) {
        char error[256];
        sprintf(error,
                "%s:(%d) key %s already exists in argument parser",
                __FILE__,
                __LINE__,
                long_form.c_str());

        throw std::range_error(error);
    }

    // add the new key
    add_arg(std::unique_ptr<arg<T>>(new arg<T>(short_form, long_form, value, help)));


    return true;
}

// setup by adding arguments
template<typename T> bool cmdargs::get_arg(const string& long_form, T& value) {
    arg<T>* argument = (arg<T>*)(*this)[long_form];

    bool is_set = argument->is_set();
    if (is_set)
        value = argument->get();


    return is_set;
}
