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
// Description: helper functions for parsers
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
#pragma once

#include <string>
#include <regex>
#include <iostream>


using namespace std;



// Parsing functions to parse string to value, on same
// interface. Functions should recognise 
template< typename T>
bool parse_data(const string & value, T & target)
{

    return false;
}

// bool template specialisation
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

// int template specialisation
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

// double template specialisation
template<>
inline bool parse_data(const string & value, double & target )
{

    regex double_regex("^((-|\\+)?[0-9]+(\\.[0-9]+)?((E|e)(-|\\+)?[0-9]+)?)$");



    std::smatch match;

    if (regex_match(value, match, double_regex))
    {
        target = std::stod(match[1]);
        //cout << "parsed [" << value << "] as [" << target << "]" << endl;

        return true;
    }

    return false;
}
    
// string template specialisation
template<>
inline bool parse_data(const string & value, string & target )
{
    target = value;
    
    return true;
}

template<typename T> inline string to_strg(T & val) {    return std::to_string(val);  }
template<> inline string to_strg(string & val) {    return val;  }
template<> inline string to_strg(bool & val) {    return val?"1":"0";  }


// helpers to check coherence of input    
template<typename T>
bool check_greater(  string name, T val, T min_)
{
    if (val > min_)
        return true;
    else
    {
        cout << "Bad range for " << name << " (cond: min("<<min_<<") < " << val << endl;
        
        return false;
    }
}


template<typename T>
bool check_range( const string & name, T val, T min_, T max_)
{
    if (val > min_ && val < max_)
        return true;
    else
    {
        cout << "Bad range for " << name << " (cond: min("<<min_<<") < "<<val<<") < max("<<max_<<"))"<<endl;
        
        return false;
    }  
}

