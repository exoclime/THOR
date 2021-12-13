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

#include <functional>
#include <tuple>

void wrap_compute_radiative_transfer(
    // prepare_compute_flux
    long dev_starflux, // in: pil
    // state variables
    long
                dev_T_lay, // out: it, pil, io, mmm, kil   (interpolated from T_int and then used as input to other funcs)
    long        dev_T_int, // in: it, pii, ioi, mmmi, kii
    long        dev_p_lay, // in: io, mmm, kil
    long        dev_p_int, // in: ioi, mmmi, kii
    const bool& interp_and_calc_flux_step,
    // direct_beam_flux
    long z_lay,
    // spectral flux loop
    bool single_walk,
    // populate_spectral_flux_noniso
    long   F_down_wg,
    long   F_up_wg,
    long   Fc_down_wg,
    long   Fc_up_wg,
    long   F_dir_wg,
    long   Fc_dir_wg,
    double delta_tau_limit,
    // integrate_flux
    long F_down_tot,
    long F_up_tot,
    long F_dir_tot,
    long F_net,
    long F_down_band,
    long F_up_band,
    long F_dir_band,
    long F_up_TOA_spectrum,
    long zenith_angle,
    long surface_albedo,
    bool surface);

void init_alfrodull();
void init_parameters(const int&    nlayer_,
                     const bool&   iso_,
                     const double& Tstar_,
                     const bool&   real_star,
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
                     const bool&   geom_zenith_corr,
                     const double& f_factor,
                     const double& w_0_limit,
                     const double& i2s_transition,
                     const bool&   debug);

void deinit_alfrodull();

void set_clouds_data(const bool& clouds_,
                     const long& cloud_opac_lay_,
                     const long& cloud_opac_int_,
                     const long& cloud_scat_cross_lay_,
                     const long& cloud_scat_cross_int_,
                     const long& g_0_tot_lay_,
                     const long& g_0_tot_int_);

void set_z_calc_function(std::function<void()>& func);

// TODO: this shouldn't be visible externally
void allocate();

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

void prepare_planck_table();
void correct_incident_energy(long starflux_array_ptr, bool real_star, bool energy_budge_correction);
