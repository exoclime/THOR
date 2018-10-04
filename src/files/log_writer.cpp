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

#include "log_writer.h"

#include <stdio.h>
#include "directories.h"
#include <iomanip>
#include <stdexcept>
#include <memory>

log_writer::log_writer(const std::string & sim_id_,
                       const std::string & output_dir_ ) :
    simulation_ID(sim_id_),
    output_dir(output_dir_)
{

}


// Check if output log file exists in output_directory and has data
// returns true if exists and found names of output files, false if not
// if exists, returns name of last written output file, last file number, and iteration number
// throws an exception if it couldn't parse the file or open the file
bool log_writer::CheckOutputLog(int & file_number, int & iteration_number, std::string & last_file)
{
    path o(output_dir);

    o /= ("esp_write_log_" + simulation_ID + ".txt");

    printf("Testing output log for file %s.\n", o.to_string().c_str());

    bool exists = path_exists(o.to_string());
    
    if (exists)
    {
        std::string line;
        std::ifstream inputfile (o.to_string().c_str());
        
        if (inputfile.is_open())
        {
            file_number = 0;
            last_file = "";
            iteration_number = 0;

            bool read_filename = false;
            
            while ( std::getline (inputfile,line) )
            {
                int line_size = line.size();
                std::unique_ptr<char[]>  scan_buf = std::unique_ptr<char[]>(new char[line_size],
                                                                              std::default_delete<char[]>());
                
                if (sscanf(line.c_str(), "#%s\n", scan_buf.get()) == 1)
                {
                    // comment, ignore
                }
                else if (sscanf(line.c_str(),"%d %d %s\n", &iteration_number, &file_number, scan_buf.get()) == 3)
                {                    
                    // data line, ok
                    read_filename = true;
                    last_file = string(scan_buf.get());
                }
                else
                {   
                    // Should not happen
                    throw std::runtime_error("Exception: unrecognised line in input file: ["+line+"].");
                }
            }


            inputfile.close();

            
            printf("Found last file iteration %d, file number %d, filename: %s\n",
                   iteration_number,
                   file_number,
                   last_file.c_str());
            

            return read_filename;
        }
        else
        {
            throw std::runtime_error("Exception: error opening input file: "+o.to_string()+".");
        }
    }
    else
    {
        printf("No log output file found.");
        
        return false;
    }
    
}


void log_writer::OpenOutputLogForWrite(bool append)
{
    path o(output_dir);

    
    o /= ("esp_write_log_" + simulation_ID + ".txt");

    printf(" Output file for result write log: %s.\n", o.to_string().c_str());
    
    if (append)
        fileoutput_output_file.open(o.to_string(),std::ofstream::out | std::ofstream::app);
    else
    {
        // start new file
        fileoutput_output_file.open(o.to_string(),std::ofstream::out);
        // append header
        fileoutput_output_file << "#"
                               << "step_number\t"
                               << "file_number\t"
                               << "filename" << std::endl;
    }
    
}

void log_writer::WriteOutputLog(int step_number, int file_number, string filename)
{
    fileoutput_output_file << step_number << "\t"
                           << file_number << "\t"
                           << filename << std::endl;

    fileoutput_output_file.flush();
}

// *************************************************************************************************
// * conservation file
int log_writer::PrepareConservationFile(bool append)
{
    path o(output_dir);

    o /= ("esp_global_" + simulation_ID + ".txt");
    //printf("Output conservation file to %s\n", o.to_string().c_str());

    // Open for write
    if (append)
    {
        // write and append, open at end of file
        conservation_output_file.open(o.to_string(), std::ofstream::out | std::ofstream::app  );
    }
    else
    {
        // start a new file
        conservation_output_file.open(o.to_string(), std::ofstream::out );

        // output file header
        conservation_output_file << std::setprecision(16);
        conservation_output_file << "#" 
                                 << "current_step" << "\t"
                                 << "simulation_time" << "\t"
                                 << "GlobalE_h" << "\t"
                                 << "GlobalMass_h" << "\t"
                                 << "GlobalAMx_h" << "\t"
                                 << "GlobalAMy_h" << "\t"
                                 << "GlobalAMz_h" << std::endl;
    }
    

    return 0;
    
};


void log_writer::OutputConservation(int current_step,
                                    double simulation_time,
                                    double GlobalE_h,
                                    double GlobalMass_h,
                                    double GlobalAMx_h,
                                    double GlobalAMy_h,
                                    double GlobalAMz_h)
{
    //printf("output conservation\n");
    
    // output global conservation values
    conservation_output_file << current_step << "\t"
                             << simulation_time << "\t"
                             << GlobalE_h << "\t"
                             << GlobalMass_h << "\t"
                             << GlobalAMx_h << "\t"
                             << GlobalAMy_h << "\t"
                             << GlobalAMz_h << std::endl;
    

    // flush file to disk
    conservation_output_file.flush();    
}

// *************************************************************************************************
// * diagnostics file
int log_writer::PrepareDiagnosticsFile(bool append)
{
    path o(output_dir);

    o /= ("esp_diagnostics_" + simulation_ID + ".txt");
    //printf("Output conservation file to %s\n", o.to_string().c_str());

    // Open for write
    if (append)
    {
        // write and append, open at end of file
        diagnostics_output_file.open(o.to_string(), std::ofstream::out | std::ofstream::app  );
    }
    else
    {
        // start a new file
        diagnostics_output_file.open(o.to_string(), std::ofstream::out );

        // output file header
        diagnostics_output_file << std::setprecision(16);
        diagnostics_output_file << "#" 
                                << "current_step" << "\t"
                                << "simulation_time" << "\t"
                                << "total_bytes" << "\t"
                                << "free_bytes" << "\t"
                                << "elapsed_time"<< "\t"
                                << "time_left" << "\t"
                                << "mean_delta_per_step" << "\t"
                                << "end_time" << std::endl;
    }
    

    return 0;
    
};


void log_writer::OutputDiagnostics(int current_step,
                                   double simulation_time,
                                   size_t total_bytes,
                                   size_t free_bytes,
                                   double elapsed_time,
                                   double time_left,
                                   double mean_delta_per_step,
                                   std::time_t end_time)
{
    //printf("output conservation\n");
    
    // output global conservation values
    diagnostics_output_file << current_step << "\t"
                             << simulation_time << "\t"
                             << total_bytes << "\t"
                             << free_bytes << "\t"
                             << elapsed_time << "\t"
                             << time_left << "\t"
                             << mean_delta_per_step << "\t"
                             << end_time << std::endl;
    

    // flush file to disk
    diagnostics_output_file.flush();    
}
