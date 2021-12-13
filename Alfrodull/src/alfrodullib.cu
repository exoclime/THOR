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
// Code to call alfrodull_engine from python and integrate into
// Helios. Was used in development to check code modifications were still working in HELIOS
// Note: not up to date, needs to be adapted to multi-column version for reuse in HELIOS
//       left here as is if there is some interest to reuse alfrodull in HELIOS
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
//
////////////////////////////////////////////////////////////////////////


#include "alfrodullib.h"

#include "integrate_flux.h"

#include "interpolate_values.h"

#include "calculate_physics.h"

#include "alfrodull_engine.h"

#include <cstdio>
#include <memory>

std::unique_ptr<alfrodull_engine> Alf_ptr = nullptr;


void init_alfrodull() {
    printf("Create Alfrodull Engine\n");

    Alf_ptr = std::make_unique<alfrodull_engine>();

    Alf_ptr->init();
}

void deinit_alfrodull() {
    cudaError_t err = cudaGetLastError();

    // Check device query
    if (err != cudaSuccess) {
        printf("deinit_alfrodull: cuda error: %s\n", cudaGetErrorString(err));
    }

    printf("Clean up Alfrodull Engine\n");
    Alf_ptr = nullptr;
}

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
                     const double& mu_star_limit,
                     const bool&   debug) {
    if (Alf_ptr == nullptr) {
        printf("ERROR: Alfrodull Engine not initialised");
        return;
    }
    const int num_cols = 1;
    Alf_ptr->set_parameters(nlayer_,
                            iso_,
                            Tstar_,
                            real_star,
                            false, // null_planck_function
                            fake_opac,
                            g_0,
                            epsi,
                            epsilon2,
                            scat,
                            scat_corr,
                            R_planet,
                            R_star,
                            a,
                            dir_beam,
                            // geom_zenith_corr,
                            f_factor,
                            w_0_limit,
                            i2s_transition,
                            mu_star_limit,
                            1,    // wiggle iteration max
                            true, // G_pm limit on full G_pm
                            num_cols,
                            debug);
}

void allocate() {
    if (Alf_ptr == nullptr) {
        printf("ERROR: Alfrodull Engine not initialised");
        return;
    }

    Alf_ptr->allocate_internal_variables();
}

// TODO: this is ugly and should not exist!
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
get_device_pointers_for_helios_write() {
    if (Alf_ptr == nullptr) {
        printf("ERROR: Alfrodull Engine not initialised");
        return std::make_tuple(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }


    return Alf_ptr->get_device_pointers_for_helios_write();
}

std::tuple<long, long, long, long, int, int> get_opac_data_for_helios() {
    if (Alf_ptr == nullptr) {
        printf("ERROR: Alfrodull Engine not initialised");
        return std::make_tuple(0, 0, 0, 0, 0, 0);
    }

    return Alf_ptr->get_opac_data_for_helios();
}


void prepare_planck_table() {
    printf("Preparing planck table\n");
    if (Alf_ptr != nullptr)
        Alf_ptr->prepare_planck_table();
    else
        printf("ERROR: prepare_planck_table: no Alf_ptr\n");
}

void correct_incident_energy(long starflux_array_ptr,
                             bool real_star,
                             bool energy_budget_correction) {
    printf("Correcting incident energy\n");

    if (Alf_ptr != nullptr)
        Alf_ptr->correct_incident_energy(
            (double*)starflux_array_ptr, real_star, energy_budget_correction);
    else
        printf("ERROR: correct_incident_energy : no Alf_ptr\n");
}


void set_z_calc_function(std::function<void()>& func) {
    if (Alf_ptr != nullptr)
        Alf_ptr->set_z_calc_func(func);
}


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
    bool surface) {
    if (Alf_ptr != nullptr) {
        int num_col = 1;
        Alf_ptr->compute_radiative_transfer(
            (double*)dev_starflux, // in: pil
            (double*)
                dev_T_lay, // out: it, pil, io, mmm, kil   (interpolated from T_int and then used as input to other funcs)
            (double*)dev_T_int, // in: it, pii, ioi, mmmi, kii
            (double*)dev_p_lay, // in: io, mmm, kil
            (double*)dev_p_int, // in: ioi, mmmi, kii
            true,
            interp_and_calc_flux_step,

            // direct_beam_flux
            (double*)z_lay,
            // spectral flux loop
            single_walk,
            // populate_spectral_flux_noniso
            (double*)F_down_wg,
            (double*)F_up_wg,
            (double*)Fc_down_wg,
            (double*)Fc_up_wg,
            (double*)F_dir_wg,
            (double*)Fc_dir_wg,
            delta_tau_limit,
            // integrate_flux
            (double*)F_down_tot,
            (double*)F_up_tot,
            (double*)F_dir_tot,
            (double*)F_net,
            (double*)F_down_band,
            (double*)F_up_band,
            (double*)F_dir_band,
            (double*)F_up_TOA_spectrum,
            (double*)zenith_angle,
            (double*)surface_albedo,
            num_col,
            1, // dummy
            surface);
    }
    else
        printf("ERROR: compute_radiative_transfer : no Alf_ptr\n");
}

void set_clouds_data(const bool&   clouds_,
                     const long&   cloud_opac_lay_,
                     const long&   cloud_opac_int_,
                     const long&   cloud_scat_cross_lay_,
                     const long&   cloud_scat_cross_int_,
                     const long&   g_0_tot_lay_,
                     const long&   g_0_tot_int_,
                     const double& fcloud_) {
    if (Alf_ptr == nullptr) {
        printf("ERROR: Alfrodull Engine not initialised");
        return;
    }

    Alf_ptr->set_clouds_data(clouds_,
                             (double*)cloud_opac_lay_,
                             (double*)cloud_opac_int_,
                             (double*)cloud_scat_cross_lay_,
                             (double*)cloud_scat_cross_int_,
                             (double*)g_0_tot_lay_,
                             (double*)g_0_tot_int_,
                             fcloud_);
}
