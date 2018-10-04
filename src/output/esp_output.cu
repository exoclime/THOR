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
// Description: Writes the output.
//
//
// Method: Uses hdf5 files to write the output.
//
//
// Known limitations: None.
//
//
// Known issues: None.
//
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "esp.h"
#include "hdf5.h"
#include "storage.h"

#include "phy_modules.h"
#include "directories.h"

#include <iomanip>

#include <fstream>
#include <stdexcept>

__host__ void ESP::CopyConservationToHost()
{
        cudaMemcpy(Etotal_h   , Etotal_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(Mass_h     , Mass_d     , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(AngMomx_h  , AngMomx_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(AngMomy_h  , AngMomy_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(AngMomz_h  , AngMomz_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::CopyGlobalToHost()
{
    // Transfer global conservation values to host
        cudaMemcpy(&GlobalE_h   , GlobalE_d    , sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&GlobalMass_h, GlobalMass_d , sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&GlobalAMx_h , GlobalAMx_d  , sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&GlobalAMy_h , GlobalAMy_d  , sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&GlobalAMz_h , GlobalAMz_d  , sizeof(double), cudaMemcpyDeviceToHost);
}


__host__ void ESP::CopyToHost()
{
//
//  Description: Transfer diagnostics from the device to the host.
//
        cudaMemcpy(Rho_h      , Rho_d      , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(Wh_h       , Wh_d       , point_num * nvi * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(pressure_h , pressure_d , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(Mh_h       , Mh_d       , 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::Output(int    fidx           , // Index of output file
                          double Cp             , // Specific heat capacities [J/(Kg K)]
                          double Rd             , // Gas constant [J/(Kg K)]
                          double Omega          , // Rotation rate [s-1]
                          double Gravit         , // Gravitational acceleration [m/s2]
                          double Mmol           , // Mean molecular mass of dry air [kg]
                          double P_Ref          , // Reference surface pressure [Pa]
                          double Top_altitude   , // Top of the model's domain [m]
                          double A              ,
                          bool   conservation   ,
                          int    hstest         ,
                          bool   SpongeLayer    ){

//
//  Description: Model output.
//
    char FILE_NAME1[160];

//  GRID OUTPUT
    if(current_step == 0){
        sprintf(FILE_NAME1, "%s/esp_output_grid_%s.h5", output_dir.c_str(), simulation_ID.c_str());

        storage s(FILE_NAME1);

        s.append_table(Altitude_h,
                       nv,
                       "/Altitude",
                       "m",
                       "Altitude");

//      Altitudeh
        s.append_table(Altitudeh_h,
                       nv+1,
                       "/Altitudeh",
                       "m",
                       "Altitude at the interfaces");

//      AreasT
        s.append_table(Altitudeh_h,
                       point_num,
                       "/areasT",
                       "m^2",
                       "Main cells areas");

//      Lon-lat grid
        s.append_table(lonlat_h,
                       2*point_num,
                       "/lonlat",
                       "-",
                       "Longitudes and latitudes");

//      Number of horizontal points
        s.append_value((double)point_num,
                       "/point_num",
                       "-",
                       "Number of grid points in one level");

//      Number of vertical layers
        s.append_value((double)nv,
                       "/nv",
                       "-",
                       "Number of vertical layers");
//      point neighbours
        s.append_table(point_local_h,
                       6*point_num,
                       "/pntloc",
                       "-",
                       "Neighbours indexes");
    }

//  PLANET
    if(current_step == 0){
        sprintf(FILE_NAME1, "%s/esp_output_planet_%s.h5", output_dir.c_str(), simulation_ID.c_str());
        storage s(FILE_NAME1);

        // glevel
        s.append_value(glevel, "/glevel", "-", "Horizontal subdivision level");
        // vlevel
        s.append_value(nv,"/vlevel","-", "Vertical subdivision level");
        // spring_dynamics
        s.append_value(spring_dynamics?1.0:0.0, "/spring_dynamics","-", "Spring dynamics");
        // spring beta
        s.append_value(spring_beta, "/spring_beta", "-","Spring Beta");
        //      A
        s.append_value(A,"/A","m","Planet radius");
        //      Rd
        s.append_value(Rd,"/Rd","J/(Kg K)","Gas constant");
        //      Omega
        s.append_value(Omega,"/Omega","1/s","Rotation rate");
        //      Gravit
        s.append_value(Gravit,"/Gravit", "m/s^2","Surface gravity");
        //      Mmol
        s.append_value(Mmol,"/Mmol","kg","Mean molecular mass of dry air");
        //      P_Ref
        s.append_value( P_Ref,"/P_Ref","Pa","Reference pressure");
        //      Top_altitude
        s.append_value(Top_altitude,"/Top_altitude", "m", "Top of the model's domain");
        //      CP
        s.append_value(Cp, "/Cp", "J/(Kg K)", "Specific heat capacity");
        //      SpongeLayer option
        s.append_value(SpongeLayer?1.0:0.0, "/SpongeLayer", "-", "Using SpongeLayer?");
        // //      hstest option
        s.append_value(hstest, "/hstest", "-", "Using benchmark forcing or RT");
        if (SpongeLayer) {
            //      nlat
            s.append_value(nlat, "/nlat", "-", "number of lat rings for sponge layer");
            //      ns
            s.append_value(ns_sponge, "/ns_sponge", "-", "Bottom of sponge layer");
            //      Rv
            s.append_value(Rv_sponge, "/Rv_sponge", "1/s", "Stength of sponge layer");
        }

        phy_modules_store_init(s);
    }

//  ESP OUTPUT
    
    sprintf(FILE_NAME1, "%s/esp_output_%s_%d.h5", output_dir.c_str(), simulation_ID.c_str(), fidx);

    storage s(FILE_NAME1);
    // step index
    s.append_value(current_step,
                   "/nstep",
                   "-",
                   "Step number");

//  Simulation time
    s.append_value(simulation_time,
                   "/simulation_time",
                   "s",
                   "Simulation time");

//  Rho
    s.append_table(Rho_h,
                   nv*point_num,
                   "/Rho",
                   "kg/m^3",
                   "Density");

//  Pressure
    s.append_table(pressure_h,
                   nv*point_num,
                   "/Pressure",
                   "Pa",
                   "Pressure");

//  Mh
    s.append_table(Mh_h,
                   nv*point_num*3,
                   "/Mh",
                   "kg m/s",
                   "Horizontal Momentum");

//  Wh
     s.append_table(Wh_h,
                   nvi*point_num,
                   "/Wh",
                   "kg m/s",
                   "Vertical Momentum");

    if (conservation == true)  {
//  Etotal at each point
     s.append_table( Etotal_h,
                     nv*point_num,
                     "/Etotal",
                     "kg m^2/s^2",
                     "Total Energy");

//  Mass at each point
     s.append_table( Mass_h,
                     nv*point_num,
                     "/Mass",
                     "kg",
                     "Mass");

//  AngMomx at each point
     s.append_table( AngMomx_h,
                     nv*point_num,
                     "/AngMomx",
                     "kg m^2/s",
                     "AngMom in X");
//
// //  AngMomy at each point
     s.append_table( AngMomy_h,
                     nv*point_num,
                     "/AngMomy",
                     "kg m^2/s",
                     "AngMom in Y");
//
// //  AngMomz at each point
     s.append_table( AngMomz_h,
                     nv*point_num,
                     "/AngMomz",
                     "kg m^2/s",
                     "AngMom in Z");

//  GlobalE (total energy over entire planet)
     s.append_value( GlobalE_h,
                     "/GlobalE",
                     "kg m^2/s^2",
                     "Global Total Energy");

//  GlobalMass (total atmospheric mass over entire planet)
      s.append_value( GlobalMass_h,
                     "/GlobalMass",
                     "kg",
                     "Global Mass");

//  GlobalAMx (total angular momentum in x direction over entire planet)
      s.append_value( GlobalAMx_h,
                     "/GlobalAMx",
                     "kg m^2/s",
                     "Global AngMomX");

//  GlobalAMy (total angular momentum in y direction over entire planet)
      s.append_value( GlobalAMy_h,
                     "/GlobalAMy",
                     "kg m^2/s",
                     "Global AngMomY");

//  GlobalAMz (total angular momentum in y direction over entire planet)
      s.append_value( GlobalAMz_h,
                     "/GlobalAMz",
                     "kg m^2/s",
                     "Global AngMomZ");
      }
      phy_modules_store(s);

      char buf[256];
      
      sprintf(buf, "esp_output_%s_%d.h5", simulation_ID.c_str(), fidx);
      // Write to output f
      logwriter.WriteOutputLog(current_step, fidx, string(buf) );
}



void ESP::OutputConservation()
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


// Store path to output and prepare output files
void ESP::SetOutputParam(const std::string & sim_id_,
                         const std::string & output_dir_ )
{
    simulation_ID = sim_id_;
    output_dir = output_dir_;
}

int ESP::PrepareConservationFile()
{
    path o(output_dir);

    o /= ("esp_global_" + simulation_ID + ".txt");
    //printf("Output conservation file to %s\n", o.to_string().c_str());

    // Open for read and write.
    // TODO: will need to handle restart
    // TODO: output header at start of file
    conservation_output_file.open(o.to_string(),std::ofstream::out
                                  //  std::ofstream::in
                                  //| std::ofstream::out
                                  //  | std::ofstream::app // no append for now
        );

    conservation_output_file << std::setprecision(16);
    conservation_output_file << "#"
                             << "current_step" << "\t"
                             << "simulation_time" << "\t"
                             << "GlobalE_h" << "\t"
                             << "GlobalMass_h" << "\t"
                             << "GlobalAMx_h" << "\t"
                             << "GlobalAMy_h" << "\t"
                             << "GlobalAMz_h" << std::endl;


    return 0;

};
