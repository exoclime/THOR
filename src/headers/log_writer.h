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
// Description: Defines the main model's parameters
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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#pragma once

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>


class log_writer
{
public:
    log_writer(const std::string& sim_id_,
               const std::string& output_dir_);


    // control file log
    bool CheckOutputLog(int& file_number, int& iteration_number, std::string& last_file);
    void OpenOutputLogForWrite(bool append);
    void WriteOutputLog(int step_number, int file_number, std::string filename);

    // conservation file log
    void OutputConservation(int    current_step,
                            double simulation_time,
                            double GlobalE_h,
                            double GlobalMass_h,
                            double GlobalAMx_h,
                            double GlobalAMy_h,
                            double GlobalAMz_h);
    int  PrepareConservationFile(bool append);


    int  PrepareDiagnosticsFile(bool append);
    void OutputDiagnostics(int         current_step,
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
