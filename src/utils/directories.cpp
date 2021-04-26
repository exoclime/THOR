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

#include "directories.h"

#include <algorithm>
#include <iostream>
#include <istream>
#include <sstream>

#include <utility>
#include <vector>

#include <cstring>

#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>

using namespace std;

// #include "cuda_device_memory.h"
#include "log_writer.h"

vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string         token;
    istringstream  tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

bool path_exists(const string &path) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) {
        //log::printf("found file %s.\n", path.c_str());

        return true;
    }

    else {
        //log::printf("not found file %s.\n", path.c_str());

        // debug stat error
        /*
        int errsv = errno;
        char buf[256];

        if (strerror_r(errsv, buf, 256) == 0)
        {
            log::printf("Error reading char: %s.\n", buf);
        }
        */
        return false;
    }
}

bool create_dir(const string &dir) {
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG) == 0)
        return true;
    else
        return false;
}


bool create_output_dir(const string &output_dir) {
    // cout << "create: "<<  output_dir<<endl;
    vector<string> directories = split(output_dir, '/');

    string path = "";

    // check for absolute path
    if (output_dir[0] == '/')
        path = "/";

    for (const auto &dir : directories) {
        //cout << "dir: " << dir << endl;
        // check for empty dir, happens with absolute path
        if (dir == "")
            continue;
        // check if directory exists
        if (path != "")
            path += string("/");
        path += dir;

        if (path_exists(path)) {
            // cout << path << " exists"<<endl;

            // continue
        }
        else {
            // cout << "[" <<  path << "] does not exist, creating"<<endl;

            // create directory
            if (create_dir(path)) {
                // cout << "created " << path << " path"<<endl;
                // continue
            }
            else {
                // creation failed, failed
                log::printf("failed creating directory %s.\n", path.c_str());
                return false;
            }
        }
    }


    return true;
}


string path::suffix() {
    if (elements.size() > 0) {
        vector<string> tokens = split(elements.back(), '.');
        return tokens.back();
    }
    else
        return "";
}

vector<string> path::suffixes() {
    vector<string> out;

    if (elements.size() > 0) {
        vector<string> tokens = split(elements.back(), '.');
        for (unsigned int i = 1; i < tokens.size(); i++)
            out.push_back(tokens[i]);
    }

    return out;
}


vector<string> path::parts() {

    return elements;
}

string path::name() {
    if (elements.size() > 0)
        return elements.back();
    else
        return "";
}

string path::stem() {
    string out("");

    if (elements.size() > 0) {
        vector<string> tokens = split(elements.back(), '.');
        if (tokens.size() > 0)
            out = tokens[0];

        for (unsigned int i = 1; i < tokens.size() - 1; i++)
            out += "." + tokens[i];
    }

    return out;
}

string path::parent() {
    string p = "";
    if (is_absolute_path)
        p += "/";

    for (unsigned int i = 0; i < elements.size() - 1; i++) {
        p += elements[i];

        if (i != elements.size() - 2) // last element
            p += "/";
    }

    return p;
}

string path::to_string() {
    string p = "";
    if (is_absolute_path)
        p += "/";

    for (unsigned int i = 0; i < elements.size(); i++) {
        p += elements[i];

        if (i != elements.size() - 1) // last element
            p += "/";
    }

    return p;
}

const char *path::c_str() {
    return to_string().c_str();
}


path::path(const string &p) {
    vector<string> tokens = split(p, '/');
    if (tokens.size() > 0 && tokens[0] == "")
        is_absolute_path = true;

    for (const auto &el : tokens) {
        if (el != "")
            elements.push_back(el);
    }
}


vector<string> get_files_in_directory(const string &dir_name) {
    vector<string> files;

    DIR *          dir;
    struct dirent *ent;
    if ((dir = opendir(dir_name.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            files.push_back(dir_name + "/" + ent->d_name);
        }
        closedir(dir);
    }
    else {
        log::printf("Could not open directory %s.\n", dir_name.c_str());
    }
    return files;
}


// match filename to output pattern
// esp_output_[PLANETNAME]_[OUTPUT_NUMBER].h5
// returns planet name and output number
bool match_output_file_numbering_scheme(const string &file_path, string &basename, int &number) {
    path p(file_path);
    basename = "";
    number   = -1;

    if (p.suffix() != "h5")
        return false;

    string filename_base = p.stem();

    // split through underscores
    vector<string> split_underscores = split(filename_base, '_');

    if (split_underscores.size() >= 4) {
        try {
            number = stoi(split_underscores.back());


        } catch (...) {
            number = -1;

            return false;
        }


        if (split_underscores[0] == "esp" && split_underscores[1] == "output" && number > -1) {
            basename = split_underscores[2];

            for (unsigned int i = 3; i < split_underscores.size() - 1; i++) {
                basename += "_" + split_underscores[i];
            }

            return true;
        }
        else
            return false;
    }
    else {
        return false;
    }
}

bool find_continue_file(string &initial_conditions,
                        string &planet_filename,
                        bool    continue_sim,
                        int &   output_file_idx) {

    // build planet filename
    path   p(initial_conditions);
    int    file_number = 0;
    string basename    = "";

    string parent_path = p.parent();

    // check existence of files
    if (!path_exists(initial_conditions)) {
        log::printf("initial condition file %s not found.\n", initial_conditions.c_str());
        return false;
    }

    // Reload correct file if we are continuing from a specific file
    if (continue_sim) {
        if (!match_output_file_numbering_scheme(initial_conditions, basename, file_number)) {
            log::printf("Loading initial conditions: "
                        "Could not recognise file numbering scheme "
                        "for input %s: (found base: %s, num: %d) \n",
                        initial_conditions.c_str(),
                        basename.c_str(),
                        file_number);
            return false;
        }

        output_file_idx = file_number;

        planet_filename = p.parent() + "/esp_output_planet_" + basename + ".h5";
    }
    else {
        //test to see if initial file matches standard output naming convention
        if (match_output_file_numbering_scheme(initial_conditions, basename, file_number)) {
            //set planet file name based on standard output naming
            planet_filename = p.parent() + "/esp_output_planet_" + basename + ".h5";
        }
        else { //if not, use old naming convention for input files
            planet_filename = p.parent() + "/" + p.stem() + "_planet.h5";
        }
    }

    if (!path_exists(planet_filename)) {
        log::printf("planet_file %s not found.\n", planet_filename.c_str());
        return false;
    }

    return true;
}

bool overwrite_check(string &output_path,
                     string &simulation_ID,
                     int     output_file_idx,
                     bool    force_overwrite) {
    // Check presence of output files
    path results(output_path);

    // Get the files in the directory
    vector<string> result_files = get_files_in_directory(results.to_string());

    // match them to name pattern, get file numbers and check numbers that are greater than the
    // restart file number
    vector<pair<string, int>> matching_name_result_files;
    for (const auto &r : result_files) {
        string basename    = "";
        int    file_number = 0;
        if (match_output_file_numbering_scheme(r, basename, file_number)) {
            if (basename == simulation_ID
                && file_number > output_file_idx) //RD not sure what this do
                matching_name_result_files.emplace_back(r, file_number);
        }
    }
    // sort them by number
    std::sort(matching_name_result_files.begin(),
              matching_name_result_files.end(),
              [](const std::pair<string, int> &left, const std::pair<string, int> &right) {
                  return left.second < right.second;
              });

    if (matching_name_result_files.size() > 0) {
        if (!force_overwrite) {
            log::printf("output files already exist and would be overwritten \n"
                        "when running simulation. \n"
                        "Files found:\n");
            for (const auto &f : matching_name_result_files)
                log::printf("\t%s\n", f.first.c_str());

            log::printf(" Aborting. \n"
                        "use --overwrite to overwrite existing files.\n");

            // cuda_device_memory_manager::get_instance().deallocate();
            //
            // exit(EXIT_FAILURE);
            return false;
        }
    }
    return true;
}
