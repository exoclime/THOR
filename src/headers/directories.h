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

#include <string>
#include <vector>

using namespace std;

bool match_output_file_numbering_scheme(const string & file_path, string & basename, int & number );

bool create_output_dir(const string & output_dir);

bool path_exists(const string & path);


class path
{
public:
    path(const string & path);

    bool is_absolute() 
    {
        return is_absolute_path;
    }

    // last file extension separated by a '.'
    string suffix();
    // vector of file extensions separated by a '.' 
    vector<string> suffixes();
    // parts
    vector<string> parts();
    
    // final part of the path
    string name();
    // final part of path, without suffix
    string stem();
    // final part without numbering separated by underscore
    //   bool numbered_stem(string & stem, int & number);
//    // match stem to wildcard
//    bool match_stem(const string & stem, string & addition);
    
    
private:

    string element_name;
    vector<string> elements;

    bool is_absolute_path = false;
    
};

