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



#include "cmdargs.h"

#include <regex>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <cstdio>

cmdargs::cmdargs()
{

}

bool cmdargs::parse(int argc, char ** argv)
{
    
    bool out = true;

    // initialise parser
    parser_state = PARSE_FOR_KEY;
    parser_state_nargs = 0;
    current_arg_interface = args.end();
    
    for (int i = 0; i < argc; i++)
    {
        out &= parser_state_machine(argv[i]);
    }
    
    return out;
}


arg_interface * cmdargs::operator[](string idx)
{
    const string s = idx;

    /*
    const int n = find_me();
std::find_if(v.begin(), v.end(),
             [n](const MyClass & m) -> bool { return m.myInt == n; });
    */
    auto && it = std::find_if(args.begin(),
                              args.end(),
                              [s](const std::unique_ptr<arg_interface>  & m)
                              -> bool { return m->is_key(s); });
    
                              
    if (it != args.end())
    {
        return it->get();
    }
    else
    {
        char error[256];
        sprintf(error,
                "%s:(%d) key not found in argument list",
                __FILE__,
                __LINE__);
        
        throw std::runtime_error(error );        
    }
}


bool cmdargs::parser_state_machine(const string & str)
{
    switch(parser_state) {
    case PARSE_FOR_KEY:
        {
            regex short_key_regex("^-[A-Za-z0-9]$");
            regex long_key_regex("^--[A-Za-z0-9]+$");
            std::smatch match;
        
            if (regex_match(str, match, short_key_regex))
            {
                
            }
            else if (regex_match(str, match, short_key_regex))
            {

            }
            
            break;
        }
    case GOT_KEY:
        {
            break;
            
        }
        

    }
    
        
   return false;
}
/*
 // search predicates for keys in vectors 
bool cmdargs::long_pred(const std::unique_ptr<arg_interface> &a)
{
    return a.long_form = search_key;
}

bool cmdargs::short_pred(const std::unique_ptr<arg_interface> &a)
{
    return a.short_form = search_key;
}
*/
    
