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
//
//
// Description:
//
// Method: -
//
// Known limitations: None
//
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

#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>


class log_writer
{
public:
    log_writer(const std::string& sim_id_, const std::string& output_dir_);


    // control file log
    bool check_output_log(int& file_number, int& iteration_number, std::string& last_file);
    void open_output_log_for_write(bool append);
    void write_output_log(int step_number, int file_number, std::string filename);

    // conservation file log
    void output_conservation(int    current_step,
                             double simulation_time,
                             double GlobalE_h,
                             double GlobalMass_h,
                             double GlobalAMx_h,
                             double GlobalAMy_h,
                             double GlobalAMz_h,
                             double GlobalEnt_h);
    int  prepare_conservation_file(bool append);


    int  prepare_diagnostics_file(bool append);
    void output_diagnostics(int         current_step,
                            double      simulation_time,
                            size_t      total_bytes,
                            size_t      free_bytes,
                            double      elapsed_time,
                            double      time_left,
                            double      mean_delta_per_step,
                            std::time_t end_time);


private:
    // output variables
    std::string simulation_ID; // name of output planet
    std::string output_dir;    // output directory

    std::fstream fileoutput_output_file;
    std::fstream conservation_output_file;
    std::fstream diagnostics_output_file;
};

class log
{
public:
    static void init_logger(std::string logfile_path, bool append = false) {

        if (pFILE != nullptr) {
            fflush(pFILE);
            fclose(pFILE);
        }

        if (append)
            pFILE = fopen(logfile_path.c_str(), "a");
        else
            pFILE = fopen(logfile_path.c_str(), "w");

        if (pFILE == nullptr) {
            std::printf("Error opening log file: %s.\n", logfile_path.c_str());
            exit(-1);
        }
    }

    template<typename... Args> static void printf(Args... args) {
        // Stop GCC from complaining about format string it can't analyse
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
        if (pFILE != nullptr)
            fprintf(pFILE, args...);

        std::printf(args...);
#pragma GCC diagnostic pop
    }

    template<typename... Args> static void printf_logonly(Args... args) {
        // Stop GCC from complaining about format string it can't analyse
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
        if (pFILE != nullptr)
            fprintf(pFILE, args...);
#pragma GCC diagnostic pop
    }


    static void flush() {
        if (pFILE != nullptr)
            fflush(pFILE);
    }


private:
    log();

    ~log() {
        if (pFILE != nullptr)
            fclose(pFILE);
    }

    static FILE* pFILE;
};
