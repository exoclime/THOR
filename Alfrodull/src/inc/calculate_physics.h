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

__device__ double E_parameter(double w0, double g0, double i2s_transition);

__global__ void trans_iso(double*       trans_wg,
                          double*       delta_tau_wg,
                          double*       M_term,
                          double*       N_term,
                          double*       P_term,
                          double*       G_plus,
                          double*       G_minus,
                          double*       delta_colmass,
                          double*       opac_wg_lay,
                          double*       cloud_opac_lay,
                          double*       meanmolmass_lay,
                          double*       scat_cross_lay,
                          double*       cloud_scat_cross_lay,
                          double*       w_0,
                          double*       g0,
                          double*       g_0_tot_lay,
                          double        g_0,
                          double        epsi,
                          double        epsilon2,
                          double*       zenith_angle_cols,
                          double*       mu_star_cols,
                          double        w_0_limit,
                          bool          scat,
                          int           nbin,
                          int           ny,
                          int           nlayer,
                          int           num_cols,
                          double        fcloud,
                          bool          clouds,
                          bool          scat_corr,
                          double        mu_star_wiggle_increment,
                          bool          G_pm_limiter,
                          double        G_pm_denom_limit_for_mu_star_wiggler,
                          bool          G_pm_limit_on_full_G_pm_,
                          bool*         hit_G_pm_limit,
                          unsigned int* columns_wiggle,
                          unsigned int* columns_wiggle_request,
                          unsigned int  wiggle_request_iterator,
                          bool          debug,
                          double        i2s_transition);

__global__ void trans_noniso(double*       trans_wg_upper,
                             double*       trans_wg_lower,
                             double*       delta_tau_wg_upper,
                             double*       delta_tau_wg_lower,
                             double*       M_upper,
                             double*       M_lower,
                             double*       N_upper,
                             double*       N_lower,
                             double*       P_upper,
                             double*       P_lower,
                             double*       G_plus_upper,
                             double*       G_plus_lower,
                             double*       G_minus_upper,
                             double*       G_minus_lower,
                             double*       delta_col_upper,
                             double*       delta_col_lower,
                             double*       opac_wg_lay,
                             double*       opac_wg_int,
                             double*       cloud_opac_lay,
                             double*       cloud_opac_int,
                             double*       meanmolmass_lay,
                             double*       meanmolmass_int,
                             double*       scat_cross_lay,
                             double*       scat_cross_int,
                             double*       cloud_scat_cross_lay,
                             double*       cloud_scat_cross_int,
                             double*       w_0_upper,
                             double*       w_0_lower,
                             double*       g0_upper,
                             double*       g0_lower,
                             double*       g_0_tot_lay,
                             double*       g_0_tot_int,
                             double        g_0,
                             double        epsi,
                             double        epsilon2,
                             double*       zenith_angle_cols,
                             double*       mu_star_cols,
                             double        w_0_limit,
                             bool          scat,
                             int           nbin,
                             int           ny,
                             int           nlayer,
                             int           num_cols,
                             double        fcloud,
                             bool          clouds,
                             bool          scat_corr,
                             double        mu_star_wiggle_increment,
                             bool          G_pm_limiter,
                             double        G_pm_denom_limit_for_mu_star_wiggler,
                             bool          G_pm_limit_on_full_G_pm_,
                             bool*         hit_G_pm_limit,
                             unsigned int* columns_wiggle,
                             unsigned int* columns_wiggle_request,
                             unsigned int  wiggle_request_iterator,
                             bool          debug,
                             double        i2s_transition,
                             int           column_idx);
