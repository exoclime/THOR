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


#include <istream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include <sys/stat.h>
#include <dirent.h>


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

    for (const auto & dir : directories)
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


string path::suffix()
{
    if (elements.size() > 0)
    {
        vector<string> tokens = split(elements.back(), '.');
        return tokens.back();        
    }
    else
        return "";
}

vector<string> path::suffixes()
{
    vector<string> out;

    if (elements.size() > 0)
    {
        vector<string> tokens = split(elements.back(), '.');
        for (unsigned int i = 1; i < tokens.size(); i++)
            out.push_back(tokens[i]);
    }
    
    return out;
}


vector<string> path::parts()
{
    
    return elements;
}

string path::name()
{
    if (elements.size() > 0)
        return elements.back();
    else
        return "";    
}

string path::stem()
{
    string out("");

    if (elements.size() > 0)
    {
        vector<string> tokens = split(elements.back(), '.');
        if (tokens.size() > 0)
            out = tokens[0];
        
        for (unsigned int i = 1; i < tokens.size() - 1; i++)
            out += "." + tokens[i];
    }
    
    return out;
}

string path::parent()
{
    string p = "";
    if (is_absolute_path)
        p += "/";

    for (unsigned int i = 0; i < elements.size() - 1; i++)
    {
        p += elements[i];
        
        if(i != elements.size() - 2) // last element
            p += "/";
            
    }
    
    return p;
    
}

string path::to_string()
{
    string p = "";
    if (is_absolute_path)
        p += "/";
        
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        p += elements[i];
        
        if(i != elements.size() - 1) // last element
            p += "/";
    }
    
    return p;
    
}

const char * path::c_str()
{
    return to_string().c_str();
}


path::path(const string & p)
{    
    vector<string> tokens = split(p, '/');
    if (tokens.size() > 0 && tokens[0] == "")
        is_absolute_path = true;
    
    for (const auto & el : tokens)
    {
        if (el != "")
            elements.push_back(el);
    }
}



vector<string> get_files_in_directory(const string & dir_name)
{
    vector<string> files;

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dir_name.c_str() ) ) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            files.push_back(dir_name + "/" + ent->d_name);
        }
        closedir (dir);
    } else {
        printf("Could not open directory %s.", dir_name.c_str());
    }
    return files;
}


// match filename to output pattern
// esp_output_[PLANETNAME]_[OUTPUT_NUMBER].h5
// returns planet name and output number
bool match_output_file_numbering_scheme(const string & file_path, string & basename, int & number )
{
    path p(file_path);
    basename = "";
    number = -1;
    
    if (p.suffix() != "h5")
        return false;

    string filename_base = p.stem();

    // split through underscores
    vector<string> split_underscores = split(filename_base, '_');

    if (split_underscores.size() >= 4 )
    {
        try
        {
            number = stoi(split_underscores.back() );
            
            
        }
        catch(...)
        {
            number = -1;
            
            return false;
        }
        
        
            
        if (split_underscores[0] == "esp"
            && split_underscores[1] == "output"
            && number > -1)
        {
            basename = split_underscores[2];

            for (unsigned int i = 3; i < split_underscores.size() -1; i++)
            {
                basename += "_" + split_underscores[i];
            }

            return true;
        }
        else
            return false;        
    }
    else
    {
        return false;
    }
}

