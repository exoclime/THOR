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
// Description: Declares functions and main variables on the host and device
//
// Method: -
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
#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include <string>

#include <fstream>
#include <iostream>

#include "define.h"
#include "log_writer.h"
#include "planet.h"

#include "dyn/phy_modules_device.h"

class ESP
{

public:
    ///////////////////////////
    //  General variable
    const int    point_num;
    const int    nv;
    const int    nvi;
    const int    nl_region;
    const int    nr;
    const int    nlat;
    const int    glevel;
    const bool   spring_dynamics;
    const double spring_beta;

    ///////////////////////////
    //  Host
    int *point_local_h;
    int *maps_h;

    double *lonlat_h;

    double *Altitude_h;
    double *Altitudeh_h;

    double *nvecoa_h;
    double *nvecti_h;
    double *nvecte_h;
    double *areasT_h;
    double *areasTr_h;

    double *div_h;
    double *grad_h;

    double *func_r_h;

    double *Rho_h;
    double *pressure_h;
    double *temperature_h;
    double *Mh_h;
    double *W_h;
    double *Wh_h;

    double *Kdhz_h;
    double *Kdh4_h;

    bool check_h;

    int *  zonal_mean_tab_h;
    double Rv_sponge;
    double ns_sponge;
    double t_shrink;

    //  energy, ang momentum and mass conservation
    double *Etotal_h;     //total energy (internal+kinetic+gravit) in control volume
    double  GlobalE_h;    //total energy over entire atmosphere
    double *Mass_h;       //mass in control volume
    double  GlobalMass_h; //total mass of atmosphere
    double *AngMomx_h;
    double *AngMomy_h;
    double *AngMomz_h;
    double  GlobalAMx_h;
    double  GlobalAMy_h;
    double  GlobalAMz_h;

    ///////////////////////////
    //  Device
    int *point_local_d;
    int *maps_d;

    double *Altitude_d;
    double *Altitudeh_d;

    double *nvecoa_d;
    double *nvecti_d;
    double *nvecte_d;
    double *areasT_d;
    double *areasTr_d;

    double *lonlat_d;

    double *div_d;
    double *grad_d;

    double *func_r_d;

    double *temperature_d;
    double *Mh_d;
    double *W_d;
    double *Wh_d;

    double *h_d;
    double *hh_d;

    double *Rho_d;
    double *pressure_d;

    double *Adv_d;

    double *SlowMh_d;
    double *SlowWh_d;
    double *SlowRho_d;
    double *Slowpressure_d;

    double *Rhos_d;
    double *pressures_d;
    double *Mhs_d;
    double *Ws_d;
    double *Whs_d;

    double *Rhok_d;
    double *pressurek_d;
    double *Mhk_d;
    double *Wk_d;
    double *Whk_d;

    double *v_d;
    double *pt_d;
    double *pth_d;

    double *gtil_d;
    double *gtilh_d;

    double *Sd_d;
    double *Sp_d;

    double *Kdhz_d;
    double *Kdh4_d;

    double *DivM_d;
    double *diffpr_d;
    double *diffmh_d;
    double *diffw_d;
    double *diffrh_d;

    double *diff_d;
    double *divg_Mh_d;
    bool *  check_d;

    double *vbar_d;
    double *vbar_h;
    double *utmp;
    double *utmp_h;
    double *vtmp;
    double *wtmp;
    int *   zonal_mean_tab_d;
    int     max_count;   // max number of points in latitude rings
    double *pressureh_d; // midpoint pressure used in dry conv adj

    //  energy, ang momentum and mass conservation
    double *Etotal_d;     //total energy (internal+kinetic+gravit) in control volume
    double *GlobalE_d;    //total energy over entire atmosphere
    double *Mass_d;       //mass in control volume
    double *GlobalMass_d; //total mass of atmosphere
    double *AngMomx_d;
    double *AngMomy_d;
    double *AngMomz_d;
    double *GlobalAMx_d;
    double *GlobalAMy_d;
    double *GlobalAMz_d;

    ///////////////////////////

    //  Functions
    // Constructor, receives all grid parameters
    ESP(int *           point_local_,
        int *           maps_,
        double *        lonlat_,
        double *        Altitude_,
        double *        Altitudeh_,
        double *        nvecoa_,
        double *        nvecti_,
        double *        nvecte_,
        double *        areasT_,
        double *        areasTr_,
        double *        div_,
        double *        grad_,
        double *        func_r_,
        int             nl_region_,
        int             nr_,
        int             nv_,
        int             nvi_,
        int             glevel_,
        bool            spring_dynamics_,
        double          spring_beta_,
        int             nlat_,
        int *           zonal_mean_tab,
        double          Rv_sponge_,
        double          ns_sponge_,
        double          t_shrink_,
        int             point_num_,
        bool            conservation,
        benchmark_types core_benchmark_,
        log_writer &    logwriter_,
        int             max_count_);

    ~ESP();

    void alloc_data(bool);

    bool initial_values(bool               rest,
                        const std::string &initial_conditions_filename,
                        const bool &       continue_sim,
                        double             timestep_dyn,
                        XPlanet &          xplanet,
                        bool               sponge,
                        bool               DeepModel,
                        int                TPprof,
                        int &              nsteps,
                        double &           simulation_start_time,
                        int &              output_file_idx,
                        bool               conservation);

    void init_timestep(int    nstep,
                       double simtime,
                       double timestep_) {
        current_step    = nstep;
        simulation_time = simtime;
        timestep        = timestep_;
    };


    void Thor(const XPlanet &Planet,   // planet parameters
              bool           HyDiff,   // Turn on/off hyper-diffusion.
              bool           DivDampP, // Turn on/off divergence damping.
              bool           NonHydro, // Turn on/off non-hydrostatic.
              bool           DeepModel);         // Turn on/off deep atmosphere.

    void ProfX(const XPlanet &Planet,
               int            conv,          //
               bool           DeepModel,     // use deep model
               int            n_out,         // output step (triggers conservation calc)
               bool           sponge,        // Use sponge layer?
               bool           shrink_sponge, // Shrink sponge after some time
               bool           conservation);           // calc/output conservation quantities


    void output(int            fidx,         // Index of output file
                const XPlanet &Planet,       // planet parameters
                bool           conservation,
                bool           SpongeLayer);


    void set_output_param(const std::string &sim_id_,
                          const std::string &output_dir_);

    void conservation(const XPlanet &planet,
                      bool           DeepModel);

    void copy_to_host();
    void copy_conservation_to_host();
    void copy_global_to_host();

    bool DeepModel;

private:
    // store if we run benchmarks
    benchmark_types core_benchmark;

    // run physics modules
    bool phy_modules_execute;


    // step counter for logging
    int    current_step;
    double simulation_time;
    double timestep;

    XPlanet planet;

    // output variables
    std::string simulation_ID; // name of output planet
    std::string output_dir;


    log_writer &logwriter;

    device_RK_array_manager phy_modules_core_arrays;
};
