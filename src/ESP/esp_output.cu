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

#include "esp.h"
#include "hdf5.h"
#include "storage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "directories.h"
#include "insolation.h"
#include "phy_modules.h"

#include <iomanip>

#include <fstream>
#include <stdexcept>

__host__ void ESP::copy_globdiag_to_host() {
    cudaMemcpy(Etotal_h, Etotal_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Entropy_h, Entropy_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Mass_h, Mass_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(AngMomx_h, AngMomx_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(AngMomy_h, AngMomy_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(AngMomz_h, AngMomz_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::copy_global_to_host() {
    // Transfer global globdiag values to host
    cudaMemcpy(&GlobalE_h, GlobalE_d, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&GlobalEnt_h, GlobalEnt_d, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&GlobalMass_h, GlobalMass_d, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&GlobalAMx_h, GlobalAMx_d, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&GlobalAMy_h, GlobalAMy_d, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&GlobalAMz_h, GlobalAMz_d, sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::copy_to_host() {
    //
    //  Description: Transfer diagnostics from the device to the host.
    //
    cudaMemcpy(Rho_h, Rho_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Wh_h, Wh_d, point_num * nvi * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(pressure_h, pressure_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Mh_h, Mh_d, 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Rd_h, Rd_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Cp_h, Cp_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(
        profx_Qheat_h, profx_Qheat_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::copy_mean_to_host() {
    //
    //  Description: Transfer mean of diagnostics from the device to the host.
    //
    cudaMemcpy(Rho_mean_h, Rho_mean_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Wh_mean_h, Wh_mean_d, point_num * nvi * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(
        pressure_mean_h, pressure_mean_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Mh_mean_h, Mh_mean_d, 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
}

__host__ void ESP::output(int                    fidx, // Index of output file
                          const SimulationSetup& sim) {

    //
    //  Description: Model output.
    //
    char FILE_NAME1[512];

    //  GRID OUTPUT
    if (current_step == 0) {
        sprintf(FILE_NAME1, "%s/esp_output_grid_%s.h5", output_dir.c_str(), simulation_ID.c_str());

        storage s(FILE_NAME1);

        s.append_table(Altitude_h, nv, "/Altitude", "m", "Altitude");

        //      Altitudeh
        s.append_table(Altitudeh_h, nv + 1, "/Altitudeh", "m", "Altitude at the interfaces");

        //      AreasT
        s.append_table(areasT_h, point_num, "/areasT", "m^2", "Main cells areas");

        //      Lon-lat grid
        s.append_table(lonlat_h, 2 * point_num, "/lonlat", "-", "Longitudes and latitudes");

        //      Number of horizontal points
        s.append_value((double)point_num, "/point_num", "-", "Number of grid points in one level");

        //      Number of vertical layers
        s.append_value((double)nv, "/nv", "-", "Number of vertical layers");

        //      point neighbours
        s.append_table(point_local_h, 6 * point_num, "/pntloc", "-", "Neighbours indexes");

        //      gradient operator
        s.append_table(grad_h, 7 * 3 * point_num, "/grad", "m^-1", "Horizontal gradient operator");

        //      divergence operator
        s.append_table(div_h, 7 * 3 * point_num, "/div", "m^-1", "Horizontal divergence operator");

        //      curl z operator
        s.append_table(
            curlz_h, 7 * 3 * point_num, "/curlz", "m^-1", "Vertical component of curl operator");
    }

    //  PLANET
    if (current_step == 0) {
        sprintf(
            FILE_NAME1, "%s/esp_output_planet_%s.h5", output_dir.c_str(), simulation_ID.c_str());
        storage s(FILE_NAME1);

        // glevel
        s.append_value(glevel, "/glevel", "-", "Horizontal subdivision level");
        // vlevel
        s.append_value(nv, "/vlevel", "-", "Vertical subdivision level");
        // spring_dynamics
        s.append_value(spring_dynamics ? 1.0 : 0.0, "/spring_dynamics", "-", "Spring dynamics");
        // spring beta
        s.append_value(spring_beta, "/spring_beta", "-", "Spring Beta");
        //      A
        s.append_value(sim.A, "/A", "m", "Planet radius");
        //      Rd
        s.append_value(sim.Rd, "/Rd", "J/(Kg K)", "Gas constant");
        //      Omega
        s.append_value(sim.Omega, "/Omega", "1/s", "Rotation rate");
        //      Gravit
        s.append_value(sim.Gravit, "/Gravit", "m/s^2", "Surface gravity");
        //      P_Ref
        s.append_value(sim.P_Ref, "/P_Ref", "Pa", "Reference pressure");
        //      Top_altitude
        s.append_value(sim.Top_altitude, "/Top_altitude", "m", "Top of the model's domain");
        //      CP
        s.append_value(sim.Cp, "/Cp", "J/(Kg K)", "Specific heat capacity");
        //      SpongeLayer option
        s.append_value(
            sim.RayleighSponge ? 1.0 : 0.0, "/RayleighSponge", "-", "Using Rayleigh SpongeLayer?");
        s.append_value(
            sim.RayleighSpongeT ? 1.0 : 0.0, "/RayleighSpongeT", "-", "Using thermal SpongeLayer?");

        s.append_value(
            sim.DiffSponge ? 1.0 : 0.0, "/DiffSponge", "-", "Using Diffusive SpongeLayer?");

        //      DeepModel option
        s.append_value(sim.DeepModel ? 1.0 : 0.0, "/DeepModel", "-", "Using Deep Model");

        //      NonHydro option
        s.append_value(
            sim.NonHydro ? 1.0 : 0.0, "/NonHydro", "-", "Using Non Hydrostatic parameter");

        //      output_mean option
        s.append_value(
            sim.output_mean ? 1.0 : 0.0, "/output_mean", "-", "outputting mean quantities");

        //      DivDampP option
        s.append_value(sim.DivDampP ? 1.0 : 0.0, "/DivDampP", "-", "Using Divergence-damping");

        //      HyDiff option
        s.append_value(sim.HyDiff ? 1.0 : 0.0, "/HyDiff", "-", "Using Hyper diffusion");

        //      Hyperdiffusion strength
        s.append_value(sim.Diffc, "/Diffc", "-", "Hyperdiffusion strength");
        s.append_value(sim.Diffc_v, "/Diffc_v", "-", "Vertical hyperdiffusion strength");

        //      Tmean
        s.append_value(sim.Tmean, "/Tmean", "K", "Mean atmospheric temperature");

        //      conv adj option
        s.append_value(sim.conv_adj ? 1.0 : 0.0, "/conv_adj", "-", "Using convection adjustment");

        //      GCM on/off option
        s.append_value(sim.gcm_off ? 1.0 : 0.0, "/gcm_off", "-", "GCM off");

        //      rest option
        s.append_value(sim.rest ? 1.0 : 0.0, "/rest", "-", "Starting from rest");

        //      core_benchmark  option
        s.append_value(
            int(core_benchmark), "/core_benchmark", "-", "Using benchmark forcing or RT");
        if (sim.RayleighSponge) {
            //      nlat
            s.append_value(
                nlat_bins, "/nlat_bins", "-", "number of lat rings for Rayleigh sponge layer");
            //      ns
            s.append_value(ns_ray_sponge, "/ns_ray_sponge", "-", "Bottom of rayleigh sponge layer");
            //      Rv
            s.append_value(Ruv_sponge,
                           "/Ruv_sponge",
                           "1/s",
                           "Strength of Rayleigh sponge layer (uv components)");
            //      Rv
            s.append_value(
                Rw_sponge, "/Rw_sponge", "1/s", "Strength of Rayleigh sponge layer (w component)");
        }
        // diff_sponge
        if (sim.DiffSponge) {
            //      order
            s.append_value(
                order_diff_sponge, "order_diff_sponge", "-", "Order of diffusive sponge");
            //      ns
            s.append_value(
                ns_diff_sponge, "/ns_diff_sponge", "-", "Bottom of diffusive sponge layer");
            //      Rv
            s.append_value(Dv_sponge, "/Dv_sponge", "-", "Strength of diffusive sponge layer");
        }

        s.append_value(Tint, "/Tint", "K", "Temperature of interior heat flux");
        s.append_value(f_lw, "/f_lw", "-", "fraction of taulw in well-mixed absorber");
        s.append_value(kappa_lw, "/kappa_lw", "-", "longwave opacity");
        s.append_value(kappa_sw, "/kappa_sw", "-", "shortwave opacity");

        // store module name in the description
        s.append_value(0.0, "/phy_module", "-", phy_modules_get_name());

        if (phy_modules_execute) {
            phy_modules_store_init(s);

            insolation.store_init(s);
        }
    }

    //  ESP OUTPUT

    sprintf(FILE_NAME1, "%s/esp_output_%s_%d.h5", output_dir.c_str(), simulation_ID.c_str(), fidx);

    storage s(FILE_NAME1);
    // step index
    s.append_value(current_step, "/nstep", "-", "Step number");

    //  Simulation time
    s.append_value(simulation_time, "/simulation_time", "s", "Simulation time");

    //  Rho
    s.append_table(Rho_h, nv * point_num, "/Rho", "kg/m^3", "Density");

    //  Pressure
    s.append_table(pressure_h, nv * point_num, "/Pressure", "Pa", "Pressure");

    //  Mh
    s.append_table(Mh_h, nv * point_num * 3, "/Mh", "kg m/s", "Horizontal Momentum");

    //  Wh
    s.append_table(Wh_h, nvi * point_num, "/Wh", "kg m/s", "Vertical Momentum");

    if (sim.globdiag == true) {
        //  Etotal at each point
        s.append_table(Etotal_h, nv * point_num, "/Etotal", "kg m^2/s^2", "Total Energy");

        s.append_table(Entropy_h, nv * point_num, "/Entropy", "kg m^2/s^2 K^-1", "Entropy");

        //  Mass at each point
        s.append_table(Mass_h, nv * point_num, "/Mass", "kg", "Mass");

        //  AngMomx at each point
        s.append_table(AngMomx_h, nv * point_num, "/AngMomx", "kg m^2/s", "AngMom in X");
        //
        // //  AngMomy at each point
        s.append_table(AngMomy_h, nv * point_num, "/AngMomy", "kg m^2/s", "AngMom in Y");
        //
        // //  AngMomz at each point
        s.append_table(AngMomz_h, nv * point_num, "/AngMomz", "kg m^2/s", "AngMom in Z");

        //  GlobalE (total energy over entire planet)
        s.append_value(GlobalE_h, "/GlobalE", "kg m^2/s^2", "Global Total Energy");

        s.append_value(GlobalEnt_h, "/GlobalEnt", "kg m^2/s^2 K^-1", "Global Entropy");

        //  GlobalMass (total atmospheric mass over entire planet)
        s.append_value(GlobalMass_h, "/GlobalMass", "kg", "Global Mass");

        //  GlobalAMx (total angular momentum in x direction over entire planet)
        s.append_value(GlobalAMx_h, "/GlobalAMx", "kg m^2/s", "Global AngMomX");

        //  GlobalAMy (total angular momentum in y direction over entire planet)
        s.append_value(GlobalAMy_h, "/GlobalAMy", "kg m^2/s", "Global AngMomY");

        //  GlobalAMz (total angular momentum in y direction over entire planet)
        s.append_value(GlobalAMz_h, "/GlobalAMz", "kg m^2/s", "Global AngMomZ");
    }

    // profX Qheat from physics modules
    s.append_table(profx_Qheat_h, nv * point_num, "/Qheat", "W m^-3", "Physics module Qheat");

    if (sim.output_mean == true) {
        //  Rho
        s.append_table(Rho_mean_h, nv * point_num, "/Rho_mean", "kg/m^3", "Mean Density");

        //  Pressure
        s.append_table(pressure_mean_h, nv * point_num, "/Pressure_mean", "Pa", "Mean Pressure");

        //  Mh
        s.append_table(
            Mh_mean_h, nv * point_num * 3, "/Mh_mean", "kg m/s", "Mean Horizontal Momentum");

        //  Wh
        s.append_table(Wh_mean_h, nvi * point_num, "/Wh_mean", "kg m/s", "Mean Vertical Momentum");
    }

    s.append_table(Rd_h, nv * point_num, "/Rd", "J/K/kg", "Local gas constant");
    s.append_table(Cp_h, nv * point_num, "/Cp", "J/K/kg", "Local heat capacity");


    if (phy_modules_execute) {
        phy_modules_store(*this, s);
        insolation.store(*this, s);
    }

    char buf[512];

    sprintf(buf, "esp_output_%s_%d.h5", simulation_ID.c_str(), fidx);
    // Write to output f
    logwriter.write_output_log(current_step, fidx, string(buf));
}

// Store path to output and prepare output files
void ESP::set_output_param(const std::string& sim_id_, const std::string& output_dir_) {
    simulation_ID = sim_id_;
    output_dir    = output_dir_;
}
