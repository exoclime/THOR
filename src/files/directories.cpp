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
// directories - helper to check and create directories
//
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//   
//
// If you use this code please cite the following reference: 
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016  
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//                     Russell Deitrick, russell.deitrick@csh.unibe.ch
//                     Urs Schroffenegger, urs.schroffenegger@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "directories.h"

#include <vector>
#include <istream>
#include <iostream>
#include <sstream>

#include <sys/stat.h>

using namespace std;



vector<string> split(const string& s, char delimiter)
{
   vector<string> tokens;
   string token;
   istringstream tokenStream(s);
   while (getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   
   return tokens;
}

bool path_exists(const string & path)
{
    struct stat st;
    if(stat(path.c_str(),&st) == 0)
        return true;
    else
        return false;
    

}

bool create_dir(const string & dir)
{
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG) == 0)
        return true;
    else
        return false;
    

}

    
    
bool create_output_dir(const string & output_dir)
{
    vector<string> directories = split(output_dir, '/');

    string path = "";

    for (auto dir : directories)
    {
        // check if directory exists
        if (path!= "")
            path += string("/");
        path += dir;

        if (path_exists(path))
        {
            //cout << path << " exists"<<endl;
            
            // continue
        }
        else
        {
            //cout << path << " does not exist, creating"<<endl;
            
            // create directory
            if (create_dir(path))
            {
                //cout << "created " << path << " path"<<endl;
                // continue
            }
            else
            {
                // creation failed, failed
                std::cout << "failed creating directory "<<path<<endl;
                return false;
                
            }
        }
        
            
            
        
    }
    
    
    return true;
}
