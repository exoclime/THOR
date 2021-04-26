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
#include <sstream>


#include "define.h"
#include "log_writer.h"
#include "simulation_setup.h"

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
    double *areas_h;

    double *div_h;
    double *grad_h;
    double *curlz_h;

    double *func_r_h;

    double *Rho_h;
    double *pressure_h;
    double *temperature_h;
    double *Mh_h;
    double *W_h;
    double *Wh_h;

    // average quantities over output interval
    double *Rho_mean_h;
    double *pressure_mean_h;
    double *Mh_mean_h;
    double *Wh_mean_h;

    //
    // double *Kdhz_h;
    // double *Kdh4_h;
    // double *Kdvz_h;
    // double *Kdv6_h;

    double *flux_vec;
    double *boundary_flux_h;
    double *boundary_flux_d;


    // guillot profile set-up, also to be borrowed by double gray scheme
    double Tint;     // temperature of internal heat flux
    double kappa_lw; //long wave opacity
    double kappa_sw; //short wave opacity
    double f_lw;     //fraction of lw optical depth due to well-mixed absorbers
    double bv_freq;  //for constbv profile, brunt-vaisala frequency

    bool check_h;

    // sponge layer settings and variables
    double t_shrink;
    bool   shrink_sponge;
    // Rayleigh sponge
    int *                 zonal_mean_tab_h;
    double                Ruv_sponge;
    double                Rw_sponge;
    double                RT_sponge;
    double                ns_ray_sponge;
    const int             nlat_bins;
    bool                  damp_uv_to_mean;
    bool                  damp_w_to_mean;
    raysp_calc_mode_types raysp_calc_mode;
    // Diffusive sponge
    double    Dv_sponge;
    double    ns_diff_sponge;
    const int order_diff_sponge;

    //  energy, ang momentum and mass globdiag
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
    double *Entropy_h;   //entropy in control volume
    double  GlobalEnt_h; //entropy over entire atmosphere

    // ultra-hot jupiter quantities
    double *Rd_h; //local value of gas constant
    double *Cp_h; //local value of heat capacity
    double *GibbsT;
    double *GibbsdG;
    int     GibbsN = 61;

    // physics module Qheat, for output
    double *profx_Qheat_h;

    ///////////////////////////
    //  Device
    int *point_local_d;
    int *maps_d;

    double *Altitude_d;
    double *Altitudeh_d;

    double *nvecoa_d;
    double *nvecti_d;
    double *nvecte_d;
    double *areasT_d;  //area of control volume (hexagon/pentagon)
    double *areasTr_d; //area of (6/5) triangles in hexagon/pentagon
    double *areas_d;   //area of sub-triangles (18 per hexagon, 15 per pentagon)

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

    // average quantities over output interval
    double *Rho_mean_d;
    double *pressure_mean_d;
    double *Mh_mean_d;
    double *Wh_mean_d;

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
    // double *pt_tau_d; // value of pot. temp. at prev. iteration in inner loop (not used x_x bad idea!)
    double *Etotal_tau_d;
    double *epotential_d;
    double *epotentialh_d;
    double *ekinetic_d;
    double *ekinetich_d;

    double *gtil_d;
    double *gtilh_d;

    double *Sd_d;
    double *Sp_d;

    double *Kdhz_d;
    double *Kdh4_d;
    double *Kdvz_d;
    double *Kdv6_d;

    double *Kdh2_d;

    double *DivM_d;
    double *diffpr_d;
    double *diffmh_d;
    double *diffw_d;
    double *diffrh_d;

    double *diffprv_d;
    double *diffmv_d;
    double *diffwv_d;
    double *diffrv_d;

    double *diff_d;
    double *diffv_d1;
    double *diffv_d2;
    double *divg_Mh_d;
    bool *  check_d;

    double *profx_Qheat_d;
    double *profx_dMh_d;
    double *profx_dWh_d;
    double *profx_dW_d;


    double *vbar_d;
    double *vbar_h;
    double *utmp;
    double *utmp_h;
    double *vtmp_h;
    double *wtmp_h;
    double *vtmp;
    double *wtmp;
    double *Tbar_d;
    double *Tbar_h;
    double *Ttmp;
    double *Ttmp_h;
    int *   zonal_mean_tab_d;
    int     max_count;   // max number of points in latitude rings
    double *pressureh_d; // midpoint pressure used in dry conv adj

    //  energy, ang momentum and mass globdiag
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
    double *Entropy_d; // entropy in control volume
    double *GlobalEnt_d;

    // ultra-hot jupiter quantities
    double *Rd_d; //local value of gas constant
    double *Cp_d; //local value of heat capacity
    double *GibbsT_d;
    double *GibbsdG_d;

    ///////////////////////////

    //  Functions
    // Constructor, receives all grid parameters
    ESP(int *                 point_local_,
        int *                 maps_,
        double *              lonlat_,
        double *              Altitude_,
        double *              Altitudeh_,
        double *              nvecoa_,
        double *              nvecti_,
        double *              nvecte_,
        double *              areasT_,
        double *              areasTr_,
        double *              areas_,
        double *              div_,
        double *              grad_,
        double *              curlz_,
        double *              func_r_,
        int                   nl_region_,
        int                   nr_,
        int                   nv_,
        int                   nvi_,
        int                   glevel_,
        bool                  spring_dynamics_,
        double                spring_beta_,
        int                   nlat_bins_,
        int *                 zonal_mean_tab,
        double                Ruv_sponge_,
        double                Rw_sponge_,
        double                RT_sponge_,
        double                ns_ray_sponge_,
        bool                  damp_uv_to_mean_,
        bool                  damp_w_to_mean_,
        raysp_calc_mode_types raysp_calc_mode_,
        double                Dv_sponge_,
        double                ns_diff_sponge_,
        int                   order_diff_sponge_,
        double                t_shrink_,
        bool                  shrink_sponge_,
        int                   point_num_,
        bool                  globdiag,
        benchmark_types       core_benchmark_,
        log_writer &          logwriter_,
        int                   max_count_,
        bool                  output_mean,
        init_PT_profile_types init_PT_profile_,
        double                Tint_,
        double                kappa_lw_,
        double                kappa_sw_,
        double                f_lw_,
        double                bv_freq_,
        uh_thermo_types       ultrahot_thermo_,
        uh_heating_types      ultrahot_heating_,
        thermo_equation_types thermo_equation_);

    ~ESP();

    void alloc_data(bool, bool);

    bool initial_values(const std::string &initial_conditions_filename,
                        const bool &       continue_sim,
                        double             timestep_dyn,
                        SimulationSetup &  sim,
                        int &              nsteps,
                        double &           simulation_start_time,
                        int &              output_file_idx);

    void init_timestep(int nstep, double simtime, double timestep_) {
        current_step    = nstep;
        simulation_time = simtime;
        timestep        = timestep_;
    };

    bool read_in_gibbs_H(int);

    // double chi_H_equilibrium(double, double);

    void Thor(const SimulationSetup &sim);

    void ProfX(const SimulationSetup &sim,
               int                    n_out); // output step (triggers globdiag calc)


    void output(int                    fidx, // Index of output file
                const SimulationSetup &sim); // planet parameters)


    void set_output_param(const std::string &sim_id_, const std::string &output_dir_);

    void globdiag(const SimulationSetup &sim);

    void copy_to_host();
    void copy_globdiag_to_host();
    void copy_global_to_host();
    void copy_mean_to_host();

    void update_mean_outputs(int);

    /// crash reporting utilities ///////////
    std::string index_to_location_scalar(int i, int first) {
        int idx, lev;
        lev = i % nv;
        idx = (i - lev) / nv;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tLev\tLat\tLong\tAlt" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << lev << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI << "\t" << Altitude_h[lev];
        return string_stream.str();
    }

    std::string index_to_location_scalar_mid(int i, int first) {
        int idx, lev;
        lev = i % nvi;
        idx = (i - lev) / nvi;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tLev\tLat\tLong\tAlt" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << lev << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI << "\t" << Altitudeh_h[lev];
        return string_stream.str();
    }

    std::string index_to_location_vector(int i, int first) {
        int idx, lev, xyz;
        xyz = i % 3;
        lev = ((i - xyz) / 3) % nv;
        idx = (i - 3 * lev - xyz) / nv / 3;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tLev\tXYZ\tLat\tLong\tAlt" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << lev << "\t" << xyz << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI << "\t" << Altitude_h[lev];
        return string_stream.str();
    }

    std::string dummy(int i, int first) {
        //placeholder function to fill the build definitions in binary_test.cpp
        std::ostringstream string_stream;
        string_stream << 0;
        return string_stream.str();
    }

    std::string index_to_location_1xn(int i, int first) {
        int idx = i;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

    std::string index_to_location_2xn(int i, int first) {
        int idx, ll;
        ll  = i % 2;
        idx = (i - ll) / 2;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tLat?\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << ll << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

    std::string index_to_location_3xn(int i, int first) {
        int idx, xyz;
        xyz = i % 3;
        idx = (i - xyz) / 3;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tXYZ\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << xyz << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

    std::string index_to_location_6xn(int i, int first) {
        int idx, side;
        side = i % 6;
        idx  = (i - side) / 6;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tSide\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << side << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

    std::string index_to_location_3x6xn(int i, int first) {
        int idx, xyz, side;
        xyz  = i % 3;
        side = ((i - xyz) / 3) % 6;
        idx  = (i - 3 * side - xyz) / 6 / 3;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tXYZ\tSide\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << xyz << "\t" << side << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

    std::string index_to_location_3x7xn(int i, int first) {
        int idx, xyz, side;
        xyz  = i % 3;
        side = ((i - xyz) / 3) % 7;
        idx  = (i - 3 * side - xyz) / 7 / 3;

        std::ostringstream string_stream;
        if (first == 1) {
            string_stream << "RawIdx\tGridIdx\tXYZ\tVertex\tLat\tLong" << std::endl;
        }
        string_stream << i << "\t" << idx << "\t" << xyz << "\t" << side << "\t"
                      << lonlat_h[idx * 2 + 1] * 180 / M_PI << "\t"
                      << lonlat_h[idx * 2] * 180 / M_PI;
        return string_stream.str();
    }

private:
    // store if we run benchmarks
    benchmark_types core_benchmark;

    init_PT_profile_types init_PT_profile;

    uh_thermo_types  ultrahot_thermo;
    uh_heating_types ultrahot_heating;

    thermo_equation_types thermo_equation;

    // run physics modules
    bool phy_modules_execute;


    // step counter for logging
    int    current_step;
    double simulation_time;
    double timestep;

    // output variables
    std::string simulation_ID; // name of output planet
    std::string output_dir;


    log_writer &logwriter;

    device_RK_array_manager phy_modules_core_arrays;
};
