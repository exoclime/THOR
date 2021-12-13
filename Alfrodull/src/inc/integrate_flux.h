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
///
// calculates the integrated upwards and downwards fluxes

__global__ void integrate_flux_band(double* F_down_wg,         // in
                                    double* F_up_wg,           // in
                                    double* F_dir_wg,          // in
                                    double* F_down_band,       // out
                                    double* F_up_band,         // out
                                    double* F_dir_band,        // out
                                    double* F_up_TOA_spectrum, // out
                                    double* gauss_weight,      // in
                                    int     nbin,
                                    int     numinterfaces,
                                    int     ny);

__global__ void integrate_flux_tot(double* deltalambda, // in
                                   double* F_down_tot,  // out
                                   double* F_up_tot,    // out
                                   double* F_dir_tot,   // out
                                   double* F_net,       // out
                                   double* F_down_band, // out
                                   double* F_up_band,   // out
                                   double* F_dir_band,  // out
                                   int     nbin,
                                   int     numinterfaces);


__global__ void integrate_flux_double(double* deltalambda,
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
                                      double* gauss_weight,
                                      int     nbin,
                                      int     numinterfaces,
                                      int     ny);

__global__ void fdir_iso(double* F_dir_wg,
                         double* planckband_lay,
                         double* delta_tau_wg,
                         double* z_lay,
                         double* mu_star,
                         double  mu_star_limit,
                         double  R_planet,
                         double  R_star,
                         double  a,
                         bool    dir_beam,
                         bool    geom_zenith_corr,
                         int     ninterface,
                         int     nbin,
                         int     ny,
                         int     num_cols);

__global__ void fdir_noniso(double* F_dir_wg,
                            double* Fc_dir_wg,
                            double* planckband_lay,
                            double* delta_tau_wg_upper,
                            double* delta_tau_wg_lower,
                            double* z_lay,
                            double* mu_star,
                            double  mu_star_limit,
                            double  R_planet,
                            double  R_star,
                            double  a,
                            bool    dir_beam,
                            bool    geom_zenith_corr,
                            int     ninterface,
                            int     nbin,
                            int     ny,
                            int     num_cols);


// calculation of the spectral fluxes, isothermal case with emphasis on on-the-fly calculations
__global__ void fband_iso_notabu(double* F_down_wg,
                                 double* F_up_wg,
                                 double* F_dir_wg,
                                 double* planckband_lay,
                                 double* w_0,
                                 double* M_term,
                                 double* N_term,
                                 double* P_term,
                                 double* G_plus,
                                 double* G_minus,
                                 double* g_0_tot,
                                 double* surface_albedo,
                                 bool    singlewalk,
                                 double  Rstar,
                                 double  a,
                                 int     numinterfaces,
                                 int     nbin,
                                 double  f_factor,
                                 double* mu_star,
                                 int     ny,
                                 int     num_cols,
                                 double  epsi,
                                 bool    dir_beam,
                                 bool    clouds,
                                 bool    scat_corr,
                                 bool    debug,
                                 double  i2s_transition);


// calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
__global__ void fband_noniso_notabu(double* F_down_wg,
                                    double* F_up_wg,
                                    double* Fc_down_wg,
                                    double* Fc_up_wg,
                                    double* F_dir_wg,
                                    double* Fc_dir_wg,
                                    double* planckband_lay,
                                    double* planckband_int,
                                    double* w_0_upper,
                                    double* w_0_lower,
                                    double* delta_tau_wg_upper,
                                    double* delta_tau_wg_lower,
                                    double* M_upper,
                                    double* M_lower,
                                    double* N_upper,
                                    double* N_lower,
                                    double* P_upper,
                                    double* P_lower,
                                    double* G_plus_upper,
                                    double* G_plus_lower,
                                    double* G_minus_upper,
                                    double* G_minus_lower,
                                    double* g_0_tot_upper,
                                    double* g_0_tot_lower,
                                    double* surface_albedo,
                                    bool    singlewalk,
                                    double  Rstar,
                                    double  a,
                                    int     numinterfaces,
                                    int     nbin,
                                    double  f_factor,
                                    double* mu_star,
                                    int     ny,
                                    int     num_cols,
                                    double  epsi,
                                    double  delta_tau_limit,
                                    bool    dir_beam,
                                    bool    clouds,
                                    bool    scat_corr,
                                    bool    debug,
                                    double  i2s_transition);


// calculation of the spectral fluxes, isothermal case, using thomas algorithm
__global__ void fband_iso_thomas(double* F_down_wg,      // out
                                 double* F_up_wg,        // out
                                 double* F_dir_wg,       // in
                                 double* planckband_lay, // in
                                 double* w_0,            // in
                                 double* M_term,         // in
                                 double* N_term,         // in
                                 double* P_term,         // in
                                 double* G_plus,         // in
                                 double* G_minus,        // in
                                 double* A_buff,         // thomas worker
                                 double* B_buff,         // thomas worker
                                 double* C_buff,         // thomas worker
                                 double* D_buff,         // thomas worker
                                 double* C_prime_buff,   // thomas worker
                                 double* D_prime_buff,   // thomas worker
                                 double* X_buff,         // thomas worker
                                 double* g_0_tot,        // in (clouds)
                                 double* surface_albedo,
                                 bool    singlewalk,
                                 double  Rstar,
                                 double  a,
                                 int     numinterfaces,
                                 int     nbin,
                                 double  f_factor,
                                 double* mu_star,
                                 int     ny,
                                 int     num_cols,
                                 double  epsi,
                                 bool    dir_beam,
                                 bool    clouds,
                                 bool    scat_corr,
                                 bool    debug,
                                 double  i2s_transition);


__global__ void fband_noniso_thomas(double* F_down_wg,
                                    double* F_up_wg,
                                    double* F_dir_wg,
                                    double* Fc_dir_wg,
                                    double* planckband_lay,
                                    double* planckband_int,
                                    double* w_0_upper,
                                    double* w_0_lower,
                                    double* delta_tau_wg_upper,
                                    double* delta_tau_wg_lower,
                                    double* M_upper,
                                    double* M_lower,
                                    double* N_upper,
                                    double* N_lower,
                                    double* P_upper,
                                    double* P_lower,
                                    double* G_plus_upper,
                                    double* G_plus_lower,
                                    double* G_minus_upper,
                                    double* G_minus_lower,
                                    double* A_buff,       // thomas worker
                                    double* B_buff,       // thomas worker
                                    double* C_buff,       // thomas worker
                                    double* D_buff,       // thomas worker
                                    double* C_prime_buff, // thomas worker
                                    double* D_prime_buff, // thomas worker
                                    double* X_buff,       // thomas worker
                                    double* g_0_tot_upper,
                                    double* g_0_tot_lower,
                                    double* surface_albedo,
                                    bool    singlewalk,
                                    double  Rstar,
                                    double  a,
                                    int     numinterfaces,
                                    int     nbin,
                                    double  f_factor,
                                    double* mu_star,
                                    int     ny,
                                    int     num_cols,
                                    double  epsi,
                                    double  delta_tau_limit,
                                    bool    dir_beam,
                                    bool    clouds,
                                    bool    scat_corr,
                                    bool    debug,
                                    double  i2s_transition);
