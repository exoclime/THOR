// ==============================================================================
// This file is part of Alfrodull.
//
//     Alfrodull is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     Alfrodull is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     Alfrodull directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Method: Helios Two Stream algorithm
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// Code contributors: Urs Schroffenegger, Matej Malik
//
// History:
// Version Date       Comment
// ======= ====       =======
// 1.0     2020-07-15 First version
//
//
////////////////////////////////////////////////////////////////////////

#pragma once


#include "cloud_opacities.h"
#include "cuda_device_memory.h"
#include "opacities.h"
#include "planck_table.h"

#include <functional>

class alfrodull_engine
{
public:
    alfrodull_engine();

    void init();

    void reset();

    void load_opacities(const string& filename, bool opacity_file_is_CGS);
    void prepare_planck_table();

    void set_parameters(const int&    nlayer_,
                        const bool&   iso_,
                        const double& T_star_,
                        const bool&   real_star,
                        const bool&   null_planck_function,
                        const double& fake_opac,
                        const double& g_0,
                        const double& epsi,
                        const double& epsilon2,
                        const bool&   scat,
                        const bool&   scat_corr,
                        const double& R_planet,
                        const double& R_star,
                        const double& a,
                        const bool&   dir_beam,
                        //const bool&   geom_zenith_corr,
                        const double& f_factor,
                        const double& w_0_limit,
                        const double& i2s_transition,
                        const double& mu_star_limit,
                        const int&    wiggle_iteration_max_,
                        const bool&   G_pm_limit_on_full_G_pm_,
                        const int&    num_parallel_columns,
                        const bool&   debug);

    void set_experimental_opacity_offset(double opac) {
        opacities.experimental_opacities_offset = opac;
    };

    bool   real_star     = false;
    double fake_opac     = 0.0;
    double g_0           = 0.0;
    double epsi          = 0.5;
    double epsilon2      = 2.0 / 3.0;
    double mu_star_limit = 0.02; // about 2 degree off
    bool   scat          = false;
    bool   scat_corr     = false;
    double R_planet      = 0.0;
    double R_star        = 0.0;
    double a             = 0.0;
    bool   dir_beam      = false;

    bool null_planck_function = false;
    bool geom_zenith_corr     = true;

    double f_factor  = 0.0;
    double w_0_limit = 0.0;

    int debug_nstep   = 0;
    int debug_col_idx = 0;

    double i2s_transition = 0.0;
    bool   debug          = false;

    bool G_pm_limiter = true;

    bool   G_pm_limit_on_full_G_pm  = true;
    double G_pm_denom_limit         = 300.0;
    double mu_star_wiggle_increment = 0.5;
    int    wiggle_iteration_max     = 10;

    int max_num_parallel_columns = 1;

    cuda_device_memory<bool> hit_G_pm_denom_limit;
    // call if using clouds, to set data array pointers
    void set_clouds_data(const bool& clouds,
                         double*     cloud_abs_cross_lay,
                         double*     cloud_abs_cross_int,
                         double*     cloud_scat_cross_lay,
                         double*     cloud_scat_cross_int,
                         double*     g_0_cloud_lay,
                         double*     g_0_cloud_int,
                         double      fcloud);

    double* cloud_abs_cross_lay  = nullptr;
    double* cloud_abs_cross_int  = nullptr;
    double* cloud_scat_cross_lay = nullptr;
    double* cloud_scat_cross_int = nullptr;
    double* g_0_cloud_lay        = nullptr;
    double* g_0_cloud_int        = nullptr;

    bool   clouds = false;
    double fcloud = 0.5;

    bool thomas = false;
    void allocate_internal_variables();

    // for prototyping wrapper for HELIOS.
    std::tuple<long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               long,
               int,
               int>
    get_device_pointers_for_helios_write();

    std::tuple<long, long, long, long, int, int> get_opac_data_for_helios();

    void correct_incident_energy(double* starflux_array_ptr,
                                 bool    real_star,
                                 bool    energy_budge_correction);

    void set_z_calc_func(std::function<void()>& fun);
    void call_z_callback();

    bool get_column_integrated_g0_w0(double* g0_, double* w0_, const int& num_cols);

    //private:
    opacity_table       opacities;
    cloud_opacity_table cloud_opacities;

    planck_table plancktable;

    // general sim parameters
    //    int nbin = 0; // should come from opacity table (?)

    int    nlayer     = 0;
    int    ninterface = 0; // nlayers + 1
    bool   iso        = false;
    double T_star     = 0.0;

    std::function<void()> calc_z_func;

    // device memory
    // mu_star: computed from zenith_angle and modified by wiggle
    cuda_device_memory<double> mu_star_cols;
    // iteration checcker for mu_star check
    // used as ping pong buffer
    cuda_device_memory<unsigned int> mu_star_iteration_buffer1;
    cuda_device_memory<unsigned int> mu_star_iteration_buffer2;
    //  scattering
    cuda_device_memory<double> scatter_cross_section_lay;
    cuda_device_memory<double> scatter_cross_section_inter;

    // planck function
    cuda_device_memory<double> planckband_lay;
    cuda_device_memory<double> planckband_int;

    // delta tau, for weights. Only used internally (on device) for flux computations
    // and shared at the end for integration over wg
    // iso
    cuda_device_memory<double> delta_tau_wg;
    // noiso
    cuda_device_memory<double> delta_tau_wg_upper;
    cuda_device_memory<double> delta_tau_wg_lower;

    //
    cuda_device_memory<double> dev_T_int;
    cuda_device_memory<double> delta_col_mass;
    cuda_device_memory<double> delta_col_upper;
    cuda_device_memory<double> delta_col_lower;
    cuda_device_memory<double> meanmolmass_int;
    cuda_device_memory<double> meanmolmass_lay;
    cuda_device_memory<double> opac_wg_lay;
    cuda_device_memory<double> opac_wg_int;
    cuda_device_memory<double> trans_wg;
    cuda_device_memory<double> trans_wg_upper;
    cuda_device_memory<double> trans_wg_lower;

    cuda_device_memory<double> gauss_weights;
    // Flux computation quantities
    // computed in trans_iso/trans_noniso
    // used in populate_spectral_flux (iso/non_iso)
    // iso
    cuda_device_memory<double> M_term;
    cuda_device_memory<double> N_term;
    cuda_device_memory<double> P_term;
    cuda_device_memory<double> G_plus;
    cuda_device_memory<double> G_minus;
    cuda_device_memory<double> w0_wg;
    cuda_device_memory<double> g0_wg;

    cuda_device_memory<double> A_buff;       // thomas worker
    cuda_device_memory<double> B_buff;       // thomas worker
    cuda_device_memory<double> C_buff;       // thomas worker
    cuda_device_memory<double> D_buff;       // thomas worker
    cuda_device_memory<double> C_prime_buff; // thomas worker
    cuda_device_memory<double> D_prime_buff; // thomas worker
    cuda_device_memory<double> X_buff;       // thomas worker
    // noniso
    cuda_device_memory<double> M_upper;
    cuda_device_memory<double> M_lower;
    cuda_device_memory<double> N_upper;
    cuda_device_memory<double> N_lower;
    cuda_device_memory<double> P_upper;
    cuda_device_memory<double> P_lower;
    cuda_device_memory<double> G_plus_upper;
    cuda_device_memory<double> G_plus_lower;
    cuda_device_memory<double> G_minus_upper;
    cuda_device_memory<double> G_minus_lower;
    cuda_device_memory<double> w0_wg_upper;
    cuda_device_memory<double> w0_wg_lower;
    cuda_device_memory<double> g0_wg_upper;
    cuda_device_memory<double> g0_wg_lower;

    // for output
    cuda_device_memory<double> w0_band;
    cuda_device_memory<double> g0_band;

    void compute_radiative_transfer(double*     dev_starflux,
                                    double*     dev_T_lay,
                                    double*     dev_T_int,
                                    double*     dev_p_lay,
                                    double*     dev_p_int,
                                    const bool& interp_temp_and_pres,
                                    const bool& interp_and_calc_flux_step,
                                    double*     z_lay,
                                    bool        single_walk,
                                    double*     F_down_wg,
                                    double*     F_up_wg,
                                    double*     Fc_down_wg,
                                    double*     Fc_up_wg,
                                    double*     F_dir_wg,
                                    double*     Fc_dir_wg,
                                    double      delta_tau_limit,
                                    double*     F_down_tot,
                                    double*     F_up_tot,
                                    double*     F_dir_tot,
                                    double*     F_net,
                                    double*     F_down_band,
                                    double*     F_up_band,
                                    double*     F_dir_band,
                                    double*     F_up_TOA_spectrum,
                                    double*     zenith_angle,
                                    double*     surface_albedo,
                                    int         num_cols,
                                    int         column_idx,
                                    bool        surface);

    bool prepare_compute_flux(double*       dev_starflux,
                              double*       dev_T_lay,
                              double*       dev_T_int,
                              double*       dev_p_lay,
                              double*       dev_p_int,
                              double*       dev_opac_wg_lay,
                              double*       dev_opac_wg_int,
                              double*       dev_meanmolmass_lay,
                              double*       dev_meanmolmass_int,
                              const bool&   real_star,
                              const double& fake_opac,
                              const bool&   interp_temp_and_pres,
                              const bool&   interp_and_calc_flux_step,
                              const int&    num_col);

    void integrate_flux(double* deltalambda,
                        double* F_down_tot,
                        double* F_up_tot,
                        double* F_dir_tot,
                        double* F_net,
                        double* F_down_wg,
                        double* F_up_wg,
                        double* F_dir_wg,
                        double* F_down_band,
                        double* F_up_band,
                        double* F_dir_band,
                        double* F_up_TOA_spectrum,
                        double* gauss_weight,
                        int     num_cols);

    void calculate_transmission_iso(double* trans_wg,             // out
                                    double* delta_colmass,        // in
                                    double* opac_wg_lay,          // in
                                    double* cloud_abs_cross_lay,  // in
                                    double* meanmolmass_lay,      // in
                                    double* cloud_scat_cross_lay, // in
                                    double* g_0_tot_lay,          // in
                                    double  g_0,
                                    double  epsi,
                                    double  epsilon_2_,
                                    double* zenith_angle_cols,
                                    bool    scat,
                                    bool    clouds,
                                    int     num_cols);

    void calculate_transmission_noniso(double* trans_wg_upper,
                                       double* trans_wg_lower,
                                       double* delta_col_upper,
                                       double* delta_col_lower,
                                       double* opac_wg_lay,
                                       double* opac_wg_int,
                                       double* cloud_abs_cross_lay,
                                       double* cloud_abs_cross_int,
                                       double* meanmolmass_lay,
                                       double* meanmolmass_int,
                                       double* cloud_scat_cross_lay,
                                       double* cloud_scat_cross_int,
                                       double* g_0_tot_lay,
                                       double* g_0_tot_int,
                                       double  g_0,
                                       double  epsi,
                                       double  epsilon_2_,
                                       double* zenith_angle_cols,
                                       bool    scat,
                                       bool    clouds,
                                       int     num_cols,
                                       int     column_idx);

    bool direct_beam_flux(double* F_dir_wg,
                          double* Fc_dir_wg,
                          double* z_lay,
                          double  R_planet,
                          double  R_star,
                          double  a,
                          bool    dir_beam,
                          bool    geom_zenith_corr,
                          int     num_cols);

    bool populate_spectral_flux_iso_thomas(double* F_down_wg, // out
                                           double* F_up_wg,   // out
                                           double* F_dir_wg,  // in
                                           double* g_0_tot,   // in
                                           double* surface_albedo,
                                           bool    singlewalk,
                                           double  Rstar,
                                           double  a,
                                           double  f_factor,
                                           double  epsi,
                                           double  w_0_limit,
                                           bool    dir_beam,
                                           bool    clouds,
                                           bool    surface,
                                           int     num_cols);

    bool populate_spectral_flux_iso(double* F_down_wg, // out
                                    double* F_up_wg,   // out
                                    double* F_dir_wg,  // in
                                    double* g_0_tot,   // in
                                    double* surface_albedo,
                                    bool    singlewalk,
                                    double  Rstar,
                                    double  a,
                                    double  f_factor,
                                    double  epsi,
                                    double  w_0_limit,
                                    bool    dir_beam,
                                    bool    clouds,
                                    bool    surface,
                                    int     num_cols);

    bool populate_spectral_flux_noniso(double* F_down_wg,
                                       double* F_up_wg,
                                       double* Fc_down_wg,
                                       double* Fc_up_wg,
                                       double* F_dir_wg,
                                       double* Fc_dir_wg,
                                       double* g_0_tot_upper,
                                       double* g_0_tot_lower,
                                       double* surface_albedo,
                                       bool    singlewalk,
                                       double  Rstar,
                                       double  a,
                                       double  f_factor,
                                       double  epsi,
                                       double  w_0_limit,
                                       double  delta_tau_limit,
                                       bool    dir_beam,
                                       bool    clouds,
                                       bool    surface,
                                       double* trans_wg_upper,
                                       double* trans_wg_lower,
                                       int     num_cols);

    bool populate_spectral_flux_noniso_thomas(double* F_down_wg,
                                              double* F_up_wg,
                                              double* F_dir_wg,
                                              double* Fc_dir_wg,
                                              double* g_0_tot_upper,
                                              double* g_0_tot_lower,
                                              double* surface_albedo,
                                              bool    singlewalk,
                                              double  Rstar,
                                              double  a,
                                              double  f_factor,
                                              double  epsi,
                                              double  w_0_limit,
                                              double  delta_tau_limit,
                                              bool    dir_beam,
                                              bool    clouds,
                                              bool    surface,
                                              double* trans_wg_upper,
                                              double* trans_wg_lower,
                                              int     num_cols);
};
