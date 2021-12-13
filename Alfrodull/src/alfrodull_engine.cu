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
// Master engine to compute physical quantities through interpolation, integrate flux
// Works on batches of comlumns
//
//
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


#include "alfrodull_engine.h"
#include "gauss_legendre_weights.h"

#include "calculate_physics.h"
#include "inc/cloud_opacities.h"
#include "integrate_flux.h"
#include "interpolate_values.h"

#include "binary_test.h"
#include "debug.h"

#include <functional>
#include <map>

#include "math_helpers.h"

using std::string;

// enable this with BENCHMARKING to have detailed benchmarking of column blocks
//#define DETAILED_BENCHMARK

#ifdef BENCHMARKING
#    ifdef DETAILED_BENCHMARK
#        define PER_BLOCK_COLUMN_BENCHMARK
#    else
#        undef USE_BENCHMARK
#        define USE_BENCHMARK(...)
#        undef BENCH_POINT_I_S
#        define BENCH_POINT_I_S(...)
#    endif
#endif // BEWNCHMARKING

void cuda_check_status_or_exit(const char* filename, int line) {
    cudaError_t err = cudaGetLastError();

    // Check device query
    if (err != cudaSuccess) {
        printf("[%s:%d] CUDA error check reports error: %s\n",
               filename,
               line,
               cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

alfrodull_engine::alfrodull_engine() {
    printf("Creating Alfrodull engine\n");
}

void alfrodull_engine::load_opacities(const string& filename, bool opacity_file_is_CGS) {
    printf("Loading opacities from %s\n", filename.c_str());

    opacities.load_opacity_table(filename, opacity_file_is_CGS);
}

void alfrodull_engine::init() {
    printf("Alfrodull Init\n");

    load_opacities("input/opac_sample.h5", /* opacity_file_is_CGS */ false);
}

void alfrodull_engine::set_parameters(const int&    nlayer_,
                                      const bool&   iso_,
                                      const double& T_star_,
                                      const bool&   real_star_,
                                      const bool&   null_planck_function_,
                                      const double& fake_opac_,
                                      const double& g_0_,
                                      const double& epsi_,
                                      const double& epsilon_2_,
                                      const bool&   scat_,
                                      const bool&   scat_corr_,
                                      const double& R_planet_,
                                      const double& R_star_,
                                      const double& a_,
                                      const bool&   dir_beam_,
                                      //const bool&   geom_zenith_corr_,
                                      const double& f_factor_,
                                      const double& w_0_limit_,
                                      const double& i2s_transition_,
                                      const double& mu_star_limit_,
                                      const int&    wiggle_iteration_max_,
                                      const bool&   G_pm_limit_on_full_G_pm_,
                                      const int&    max_num_parallel_columns_,
                                      const bool&   debug_) {
    nlayer     = nlayer_;
    ninterface = nlayer + 1;
    iso        = iso_;
    T_star     = T_star_;

    real_star = real_star_;
    fake_opac = fake_opac_;
    g_0       = g_0_;
    epsi      = epsi_;
    epsilon2  = epsilon_2_;
    scat      = scat_;
    scat_corr = scat_corr_;
    R_planet  = R_planet_;
    R_star    = R_star_;
    a         = a_;
    dir_beam  = dir_beam_;
    geom_zenith_corr =
        dir_beam_; // should always be true when dir_beam = true and is irrelevant when dir_beam = false
    f_factor  = f_factor_;
    w_0_limit = w_0_limit_;

    null_planck_function = null_planck_function_;

    i2s_transition           = i2s_transition_;
    debug                    = debug_;
    mu_star_limit            = mu_star_limit_;
    wiggle_iteration_max     = wiggle_iteration_max_;
    G_pm_limit_on_full_G_pm  = G_pm_limit_on_full_G_pm_;
    max_num_parallel_columns = max_num_parallel_columns_;
}

void alfrodull_engine::allocate_internal_variables() {
    int nlayer_nbin        = nlayer * opacities.nbin;
    int nlayer_plus2_nbin  = (nlayer + 2) * opacities.nbin;
    int ninterface_nbin    = ninterface * opacities.nbin;
    int nlayer_wg_nbin     = nlayer * opacities.ny * opacities.nbin;
    int ninterface_wg_nbin = ninterface * opacities.ny * opacities.nbin;
    int num_cols           = max_num_parallel_columns;
    if (iso) {
        delta_tau_wg.allocate(num_cols * nlayer_wg_nbin);
    }
    else {
        delta_tau_wg_upper.allocate(num_cols * nlayer_wg_nbin);
        delta_tau_wg_lower.allocate(num_cols * nlayer_wg_nbin);
    }

    if (thomas) {
        if (iso) {
            A_buff.allocate(num_cols * ninterface_wg_nbin * 4);       // thomas worker
            B_buff.allocate(num_cols * ninterface_wg_nbin * 4);       // thomas worker
            C_buff.allocate(num_cols * ninterface_wg_nbin * 4);       // thomas worker
            D_buff.allocate(num_cols * ninterface_wg_nbin * 2);       // thomas worker
            C_prime_buff.allocate(num_cols * ninterface_wg_nbin * 4); // thomas worker
            D_prime_buff.allocate(num_cols * ninterface_wg_nbin * 2); // thomas worker
            X_buff.allocate(num_cols * ninterface_wg_nbin * 2);       // thomas worker
        }
        else {
            int num_th_layers             = nlayer * 2;
            int num_th_interfaces         = num_th_layers + 1;
            int num_th_interfaces_wg_nbin = num_th_interfaces * opacities.ny * opacities.nbin;
            A_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 4);       // thomas worker
            B_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 4);       // thomas worker
            C_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 4);       // thomas worker
            D_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 2);       // thomas worker
            C_prime_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 4); // thomas worker
            D_prime_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 2); // thomas worker
            X_buff.allocate(num_cols * num_th_interfaces_wg_nbin * 2);       // thomas worker
        }
    }
    // flux computation internal quantities
    if (iso) {
        M_term.allocate(num_cols * nlayer_wg_nbin);
        N_term.allocate(num_cols * nlayer_wg_nbin);
        P_term.allocate(num_cols * nlayer_wg_nbin);
        G_plus.allocate(num_cols * nlayer_wg_nbin);
        G_minus.allocate(num_cols * nlayer_wg_nbin);
        w0_wg.allocate(num_cols * nlayer_wg_nbin);
        g0_wg.allocate(num_cols * nlayer_wg_nbin);

        g0_band.allocate(num_cols * nlayer_nbin);
        w0_band.allocate(num_cols * nlayer_nbin);

        delta_col_mass.allocate(num_cols * nlayer);
    }
    else {
        M_upper.allocate(num_cols * nlayer_wg_nbin);
        M_lower.allocate(num_cols * nlayer_wg_nbin);
        N_upper.allocate(num_cols * nlayer_wg_nbin);
        N_lower.allocate(num_cols * nlayer_wg_nbin);
        P_upper.allocate(num_cols * nlayer_wg_nbin);
        P_lower.allocate(num_cols * nlayer_wg_nbin);
        G_plus_upper.allocate(num_cols * nlayer_wg_nbin);
        G_plus_lower.allocate(num_cols * nlayer_wg_nbin);
        G_minus_upper.allocate(num_cols * nlayer_wg_nbin);
        G_minus_lower.allocate(num_cols * nlayer_wg_nbin);

        // for computations
        g0_wg_upper.allocate(num_cols * nlayer_wg_nbin);
        g0_wg_lower.allocate(num_cols * nlayer_wg_nbin);
        w0_wg_upper.allocate(num_cols * nlayer_wg_nbin);
        w0_wg_lower.allocate(num_cols * nlayer_wg_nbin);

        // used to store layer value for output
        w0_wg.allocate(num_cols * nlayer_wg_nbin);
        g0_wg.allocate(num_cols * nlayer_wg_nbin);

        g0_band.allocate(num_cols * nlayer_nbin);
        w0_band.allocate(num_cols * nlayer_nbin);

        delta_col_upper.allocate(num_cols * nlayer);
        delta_col_lower.allocate(num_cols * nlayer);
    }


    mu_star_iteration_buffer1.allocate(num_cols);
    mu_star_iteration_buffer2.allocate(num_cols);

    meanmolmass_lay.allocate(num_cols * nlayer);
    // scatter cross section layer and interface
    // those are shared for print out
    scatter_cross_section_lay.allocate(num_cols * nlayer_nbin);
    planckband_lay.allocate(num_cols * nlayer_plus2_nbin);
    opac_wg_lay.allocate(num_cols * nlayer_wg_nbin);

    if (!iso) {
        meanmolmass_int.allocate(num_cols * ninterface);

        scatter_cross_section_inter.allocate(num_cols * ninterface_nbin);
        planckband_int.allocate(num_cols * ninterface_nbin);
        opac_wg_int.allocate(num_cols * ninterface_wg_nbin);
    }


    if (iso) {
        trans_wg.allocate(num_cols * nlayer_wg_nbin);
    }
    else {
        trans_wg_upper.allocate(num_cols * nlayer_wg_nbin);
        trans_wg_lower.allocate(num_cols * nlayer_wg_nbin);
    }

    mu_star_cols.allocate(num_cols);
    hit_G_pm_denom_limit.allocate(1);

    std::unique_ptr<double[]> weights = std::make_unique<double[]>(100);
    for (int i = 0; i < opacities.ny; i++)
        weights[i] = gauss_legendre_weights[opacities.ny - 1][i];

    gauss_weights.allocate(opacities.ny);
    gauss_weights.put(weights);

    USE_BENCHMARK();

#ifdef PER_BLOCK_COLUMN_BENCHMARK
    if (iso) {
        std::map<string, output_def> debug_arrays = {
            {"meanmolmass_lay",
             {meanmolmass_lay.ptr_ref(),
              (int)meanmolmass_lay.get_size(),
              "meanmolmass_lay",
              "mmml",
              true,
              dummy}},
            {"planckband_lay",
             {planckband_lay.ptr_ref(),
              (int)planckband_lay.get_size(),
              "planckband_lay",
              "plkl",
              true,
              dummy}},
            {"planckband_int",
             {planckband_int.ptr_ref(),
              (int)planckband_int.get_size(),
              "planckband_int",
              "plki",
              true,
              dummy}},
            {"opac_wg_lay",
             {opac_wg_lay.ptr_ref(),
              (int)opac_wg_lay.get_size(),
              "opac_wg_lay",
              "opc",
              true,
              dummy}},
            {"trans_wg",
             {trans_wg.ptr_ref(), (int)trans_wg.get_size(), "trans_wg", "tr", true, dummy}},
            {"scat_cs_lay",
             {scatter_cross_section_lay.ptr_ref(),
              (int)scatter_cross_section_lay.get_size(),
              "scat_cs_lay",
              "scsl",
              true,
              dummy}},
            {"delta_tau_wg",
             {delta_tau_wg.ptr_ref(),
              (int)delta_tau_wg.get_size(),
              "delta_tau_wg",
              "dtw",
              true,
              dummy}},
            {"delta_colmass",
             {delta_col_mass.ptr_ref(),
              (int)delta_col_mass.get_size(),
              "delta_colmass",
              "dcm",
              true,
              dummy}},
            {"M_term", {M_term.ptr_ref(), (int)M_term.get_size(), "M_term", "Mt", true, dummy}},
            {"N_term", {N_term.ptr_ref(), (int)N_term.get_size(), "N_term", "Nt", true, dummy}},
            {"P_term", {P_term.ptr_ref(), (int)P_term.get_size(), "P_term", "Pt", true, dummy}},
            {"G_plus", {G_plus.ptr_ref(), (int)G_plus.get_size(), "G_plus", "Gp", true, dummy}},
            {"G_minus", {G_minus.ptr_ref(), (int)G_minus.get_size(), "G_minus", "Gm", true, dummy}},
            {"w_0", {w0_wg.ptr_ref(), (int)w0_wg.get_size(), "w_0", "w0", true, dummy}},
            {"g_0", {g0_wg.ptr_ref(), (int)g0_wg.get_size(), "g_0", "g0", true, dummy}},

        };

        BENCH_POINT_REGISTER_PHY_VARS(debug_arrays, (), ());
    }
    else {
        std::map<string, output_def> debug_arrays = {
            {"meanmolmass_lay",
             {meanmolmass_lay.ptr_ref(),
              (int)meanmolmass_lay.get_size(),
              "meanmolmass_lay",
              "mmml",
              true,
              dummy}},
            {"meanmolmass_int",
             {meanmolmass_int.ptr_ref(),
              (int)meanmolmass_int.get_size(),
              "meanmolmass_int",
              "mmmi",
              true,
              dummy}},
            {"planckband_lay",
             {planckband_lay.ptr_ref(),
              (int)planckband_lay.get_size(),
              "planckband_lay",
              "plkl",
              true,
              dummy}},
            {"planckband_int",
             {planckband_int.ptr_ref(),
              (int)planckband_int.get_size(),
              "planckband_int",
              "plki",
              true,
              dummy}},
            {"opac_wg_lay",
             {opac_wg_lay.ptr_ref(),
              (int)opac_wg_lay.get_size(),
              "opac_wg_lay",
              "opcl",
              true,
              dummy}},
            {"opac_wg_int",
             {opac_wg_int.ptr_ref(),
              (int)opac_wg_int.get_size(),
              "opac_wg_int",
              "opci",
              true,
              dummy}},
            {"trans_wg_upper",
             {trans_wg_upper.ptr_ref(),
              (int)trans_wg_upper.get_size(),
              "trans_wg_upper",
              "tru",
              true,
              dummy}},
            {"trans_wg_lower",
             {trans_wg_lower.ptr_ref(),
              (int)trans_wg_lower.get_size(),
              "trans_wg_lower",
              "trl",
              true,
              dummy}},
            {"scat_cs_lay",
             {scatter_cross_section_lay.ptr_ref(),
              (int)scatter_cross_section_lay.get_size(),
              "scat_cs_lay",
              "scsl",
              true,
              dummy}},
            {"scat_cs_int",
             {scatter_cross_section_inter.ptr_ref(),
              (int)scatter_cross_section_inter.get_size(),
              "scat_cs_int",
              "scsi",
              true,
              dummy}},
            {"delta_tau_wg_upper",
             {delta_tau_wg_upper.ptr_ref(),
              (int)delta_tau_wg_upper.get_size(),
              "delta_tau_wg_upper",
              "dtwu",
              true,
              dummy}},
            {"delta_tau_wg_lower",
             {delta_tau_wg_lower.ptr_ref(),
              (int)delta_tau_wg_lower.get_size(),
              "delta_tau_wg_lower",
              "dtwl",
              true,
              dummy}},

            {"delta_col_upper",
             {delta_col_upper.ptr_ref(),
              (int)delta_col_upper.get_size(),
              "delta_col_upper",
              "dcu",
              true,
              dummy}},
            {"delta_col_lower",
             {delta_col_lower.ptr_ref(),
              (int)delta_col_lower.get_size(),
              "delta_col_lower",
              "dcl",
              true,
              dummy}},
            {"M_upper", {M_upper.ptr_ref(), (int)M_upper.get_size(), "M_upper", "Mu", true, dummy}},
            {"M_lower", {M_lower.ptr_ref(), (int)M_lower.get_size(), "M_lower", "Ml", true, dummy}},
            {"N_upper", {N_upper.ptr_ref(), (int)N_upper.get_size(), "N_upper", "Nu", true, dummy}},
            {"N_lower", {N_lower.ptr_ref(), (int)N_lower.get_size(), "N_lower", "Nl", true, dummy}},
            {"P_upper", {P_upper.ptr_ref(), (int)P_upper.get_size(), "P_upper", "Pu", true, dummy}},
            {"P_lower", {P_lower.ptr_ref(), (int)P_lower.get_size(), "P_lower", "Pl", true, dummy}},
            {"G_plus_upper",
             {G_plus_upper.ptr_ref(),
              (int)G_plus_upper.get_size(),
              "G_plus_upper",
              "Gpu",
              true,
              dummy}},
            {"G_plus_lower",
             {G_plus_lower.ptr_ref(),
              (int)G_plus_lower.get_size(),
              "G_plus_lower",
              "Gpl",
              true,
              dummy}},
            {"G_minus_upper",
             {G_minus_upper.ptr_ref(),
              (int)G_minus_upper.get_size(),
              "G_minus_upper",
              "Gmu",
              true,
              dummy}},
            {"G_minus_lower",
             {G_minus_lower.ptr_ref(),
              (int)G_minus_lower.get_size(),
              "G_minus_lower",
              "Gml",
              true,
              dummy}},


            {"w_0_upper",
             {w0_wg_upper.ptr_ref(), (int)w0_wg_upper.get_size(), "w_0_upper", "w0u", true, dummy}},
            {"w_0_lower",
             {w0_wg_lower.ptr_ref(), (int)w0_wg_lower.get_size(), "w_0_lower", "w0l", true, dummy}},

        };
        BENCH_POINT_REGISTER_PHY_VARS(debug_arrays, (), ());
    }
#endif // PER_BLOCK_COLUMN_BENCHMARK
}

// set internal arrays to zero before loop
void alfrodull_engine::reset() {
    // delta tau, for weights. Only used internally (on device) for flux computations
    // and shared at the end for integration over wg
    if (iso) {
        delta_col_mass.zero();
        delta_tau_wg.zero();
        trans_wg.zero();
        M_term.zero();
        N_term.zero();
        P_term.zero();
        G_plus.zero();
        G_minus.zero();
        w0_wg.zero();
    }
    else {
        // noiso
        delta_tau_wg_upper.zero();
        delta_tau_wg_lower.zero();
        delta_col_upper.zero();
        delta_col_lower.zero();
        trans_wg_upper.zero();
        trans_wg_lower.zero();
        M_upper.zero();
        M_lower.zero();
        N_upper.zero();
        N_lower.zero();
        P_upper.zero();
        P_lower.zero();
        G_plus_upper.zero();
        G_plus_lower.zero();
        G_minus_upper.zero();
        G_minus_lower.zero();
        w0_wg_upper.zero();
        w0_wg_lower.zero();
    }

    dev_T_int.zero();

    meanmolmass_lay.zero();
    opac_wg_lay.zero();
    scatter_cross_section_lay.zero();
    planckband_lay.zero();

    if (!iso) {
        meanmolmass_int.zero();
        opac_wg_int.zero();
        scatter_cross_section_inter.zero();
        // planck function
        planckband_int.zero();
    }

    if (thomas) {
        A_buff.zero();       // thomas worker
        B_buff.zero();       // thomas worker
        C_buff.zero();       // thomas worker
        D_buff.zero();       // thomas worker
        C_prime_buff.zero(); // thomas worker
        D_prime_buff.zero(); // thomas worker
        X_buff.zero();       // thomas worker
    }
}

// return device pointers for helios data save
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
alfrodull_engine::get_device_pointers_for_helios_write() {
    return std::make_tuple((long)*scatter_cross_section_lay,
                           (long)*scatter_cross_section_inter,
                           (long)*opac_wg_lay,
                           (long)*planckband_lay,
                           (long)*planckband_int,
                           (long)*plancktable.planck_grid,
                           (long)*delta_tau_wg,
                           (long)*delta_tau_wg_upper,
                           (long)*delta_tau_wg_lower,
                           (long)*delta_col_mass,
                           (long)*delta_col_upper,
                           (long)*delta_col_lower,
                           (long)*meanmolmass_lay,
                           (long)*trans_wg,
                           (long)*trans_wg_upper,
                           (long)*trans_wg_lower,
                           plancktable.dim,
                           plancktable.step);
}

// get opacity data for helios
std::tuple<long, long, long, long, int, int> alfrodull_engine::get_opac_data_for_helios() {
    return std::make_tuple((long)*opacities.dev_opac_wave,
                           (long)*opacities.dev_opac_interwave,
                           (long)*opacities.dev_opac_deltawave,
                           (long)*opacities.dev_opac_y,
                           opacities.nbin,
                           opacities.ny);
}


void alfrodull_engine::prepare_planck_table() {
    plancktable.construct_planck_table(
        *opacities.dev_opac_interwave, *opacities.dev_opac_deltawave, opacities.nbin, T_star);
}

void alfrodull_engine::correct_incident_energy(double* starflux_array_ptr,
                                               bool    real_star,
                                               bool    energy_budget_correction) {
    printf("T_star %g, energy budget_correction: %s\n",
           T_star,
           energy_budget_correction ? "true" : "false");
    if (T_star > 10 && energy_budget_correction) {
        dim3 grid((int(opacities.nbin) + 15) / 16, 1, 1);
        dim3 block(16, 1, 1);

        corr_inc_energy<<<grid, block>>>(*plancktable.planck_grid,
                                         starflux_array_ptr,
                                         *opacities.dev_opac_deltawave,
                                         real_star,
                                         opacities.nbin,
                                         T_star,
                                         plancktable.dim);

        cudaDeviceSynchronize();
    }
}


void alfrodull_engine::set_z_calc_func(std::function<void()>& fun) {
    calc_z_func = fun;
}

void alfrodull_engine::call_z_callback() {
    if (calc_z_func)
        calc_z_func();
}

void alfrodull_engine::set_clouds_data(const bool& clouds_,
                                       double*     cloud_abs_cross_lay_,
                                       double*     cloud_abs_cross_int_,
                                       double*     cloud_scat_cross_lay_,
                                       double*     cloud_scat_cross_int_,
                                       double*     g_0_cloud_lay_,
                                       double*     g_0_cloud_int_,
                                       double      fcloud_) {
    // For now, cloud data are values per wavelength bins, with the input computed for the
    // correct wavelength bins, and used for the full column
    // so this is directly forwarded as "layer and interface" values,
    // for usage per volume, we need to allocate these arrays per volume
    // and an interpolation/lookup function from the input table to the volume
    // element parameters (P,T,...)  needs to be implemented, similar to opacity lookup.
    cloud_abs_cross_lay  = cloud_abs_cross_lay_;
    cloud_abs_cross_int  = cloud_abs_cross_int_;
    cloud_scat_cross_lay = cloud_scat_cross_lay_;
    cloud_scat_cross_int = cloud_scat_cross_int_;
    g_0_cloud_lay        = g_0_cloud_lay_;
    g_0_cloud_int        = g_0_cloud_int_;

    clouds = clouds_;
    fcloud = fcloud_;
}


// var already present:
// bool iso
void alfrodull_engine::compute_radiative_transfer(
    double* dev_starflux, // in: pil
    double*
                dev_T_lay_cols, // out: it, pil, io, mmm, kil   (interpolated from T_int and then used as input to other funcs)
    double*     dev_T_int_cols, // in: it, pii, ioi, mmmi, kii
    double*     dev_p_lay_cols, // in: io, mmm, kil
    double*     dev_p_int_cols, // in: ioi, mmmi, kii
    const bool& interpolate_temp_and_pres,
    const bool& interp_and_calc_flux_step,
    double*     z_lay,
    bool        single_walk,
    double*     F_down_wg_cols,
    double*     F_up_wg_cols,
    double*     Fc_down_wg_cols,
    double*     Fc_up_wg_cols,
    double*     F_dir_wg_cols,
    double*     Fc_dir_wg_cols,
    double      delta_tau_limit,
    double*     F_down_tot_cols,
    double*     F_up_tot_cols,
    double*     F_dir_tot_cols,
    double*     F_net_cols,
    double*     F_down_band_cols,
    double*     F_up_band_cols,
    double*     F_dir_band_cols,
    double*     F_up_TOA_spectrum_cols,
    double*     zenith_angle_cols,
    double*     surface_albedo,
    int         num_cols,
    int         column_idx,
    bool        surface) // number of columns this function works on
{
    USE_BENCHMARK();
    {
        prepare_compute_flux(
            dev_starflux,
            dev_T_lay_cols, // out: it, pil, io, mmm, kil   (interpolated from T_int and then used as input to other funcs)
            dev_T_int_cols,   // in: it, pii, ioi, mmmi, kii
            dev_p_lay_cols,   // in: io, mmm, kil
            dev_p_int_cols,   // in: ioi, mmmi, kii
            *opac_wg_lay,     // out: io
            *opac_wg_int,     // out: ioi
            *meanmolmass_lay, // out: mmm
            *meanmolmass_int, // out: mmmi
            real_star,        // pil
            fake_opac,        // io
            interpolate_temp_and_pres,
            interp_and_calc_flux_step,
            num_cols);

        cuda_check_status_or_exit(__FILE__, __LINE__);


        BENCH_POINT_I_S(debug_nstep,
                        debug_col_idx,
                        "Alf_prep_flx",
                        (),
                        ("opac_wg_lay",
                         "meanmolmass_lay",
                         "cloud_scat_cross_lay",
                         "scat_cs_lay",
                         "planckband_lay"));
    }
    double* deltalambda = *opacities.dev_opac_deltawave;


    // also lookup and interpolate cloud values here if cloud values
    // per volume element is needed
    // fill in g_0_cloud_lay, g_0_cloud_int, cloud_scat_cross_lay, cloud_scat_cross_int
    if (interp_and_calc_flux_step) {
        if (iso) {
            BENCH_POINT_I_S(debug_nstep, debug_col_idx, "Alf_prep_II", (), ("delta_colmass"));

            calculate_transmission_iso(*trans_wg,            // out
                                       *delta_col_mass,      // in
                                       *opac_wg_lay,         // in
                                       cloud_abs_cross_lay,  // in
                                       *meanmolmass_lay,     // in
                                       cloud_scat_cross_lay, // in
                                       g_0_cloud_lay,        // in
                                       g_0,
                                       epsi,
                                       epsilon2,
                                       zenith_angle_cols,
                                       scat,
                                       clouds,
                                       num_cols);

            BENCH_POINT_I_S(debug_nstep,
                            debug_col_idx,
                            "Alf_comp_trans",
                            (),
                            ("delta_colmass",
                             "trans_wg",
                             "delta_tau_wg",
                             "M_term",
                             "N_term",
                             "P_term",
                             "w_0",
                             "g_0",
                             "G_plus",
                             "G_minus",
                             "opac_wg_lay",
                             "meanmolmass_lay",
                             "scat_cs_lay", ));
        }
        else {
            BENCH_POINT_I_S(debug_nstep,
                            debug_col_idx,
                            "Alf_prep_II",
                            (),
                            ("delta_col_upper",
                             "delta_col_lower",
                             "meanmolmass_int",
                             "scat_cs_int",
                             "opac_wg_int"));
            calculate_transmission_noniso(*trans_wg_upper,
                                          *trans_wg_lower,
                                          *delta_col_upper,
                                          *delta_col_lower,
                                          *opac_wg_lay,
                                          *opac_wg_int,
                                          cloud_abs_cross_lay,
                                          cloud_abs_cross_int,
                                          *meanmolmass_lay,
                                          *meanmolmass_int,
                                          cloud_scat_cross_lay,
                                          cloud_scat_cross_int,
                                          g_0_cloud_lay,
                                          g_0_cloud_int,
                                          g_0,
                                          epsi,
                                          epsilon2,
                                          zenith_angle_cols,
                                          scat,
                                          clouds,
                                          num_cols,
                                          column_idx);
            BENCH_POINT_I_S(debug_nstep,
                            debug_col_idx,
                            "Alf_comp_trans",
                            (),
                            ("trans_wg_upper",
                             "trans_wg_lower",
                             "delta_tau_wg_upper",
                             "delta_tau_wg_lower",
                             "planckband_lay",
                             "planckband_int",
                             "M_upper",
                             "M_lower",
                             "N_upper",
                             "N_lower",
                             "P_upper",
                             "P_lower",
                             "G_plus_upper",
                             "G_plus_lower",
                             "G_minus_upper",
                             "G_minus_lower",
                             "w_0_upper",
                             "w_0_lower",
                             "opac_wg_lay",
                             "opac_wg_int",
                             "meanmolmass_lay",
                             "meanmolmass_int",
                             "scat_cs_lay",
                             "scat_cs_int"));
        }

        cuda_check_status_or_exit(__FILE__, __LINE__);
        call_z_callback();

        direct_beam_flux(F_dir_wg_cols,
                         Fc_dir_wg_cols,
                         z_lay,
                         R_planet,
                         R_star,
                         a,
                         dir_beam,
                         geom_zenith_corr,
                         num_cols);

        BENCH_POINT_I_S(debug_nstep, debug_col_idx, "Alf_dir_beam_trans", (), ("F_dir_wg"));

        cuda_check_status_or_exit(__FILE__, __LINE__);
    }

    if (thomas) {
        if (iso) {
            populate_spectral_flux_iso_thomas(F_down_wg_cols, // out
                                              F_up_wg_cols,   // out
                                              F_dir_wg_cols,  // in
                                              *g0_wg,         // in
                                              surface_albedo,
                                              single_walk,
                                              R_star,
                                              a,
                                              f_factor,
                                              epsi,
                                              w_0_limit,
                                              dir_beam,
                                              clouds,
                                              surface,
                                              num_cols);
        }
        else {
            populate_spectral_flux_noniso_thomas(F_down_wg_cols,
                                                 F_up_wg_cols,
                                                 F_dir_wg_cols,
                                                 Fc_dir_wg_cols,
                                                 *g0_wg_upper,
                                                 *g0_wg_lower,
                                                 surface_albedo,
                                                 single_walk,
                                                 R_star,
                                                 a,
                                                 f_factor,
                                                 epsi,
                                                 w_0_limit,
                                                 delta_tau_limit,
                                                 dir_beam,
                                                 clouds,
                                                 surface,
                                                 *trans_wg_upper,
                                                 *trans_wg_lower,
                                                 num_cols);
        }
        cuda_check_status_or_exit(__FILE__, __LINE__);

        BENCH_POINT_I_S(
            debug_nstep, debug_col_idx, "Alf_pop_spec_flx_thomas", (), ("F_up_wg", "F_down_wg"));
    }
    else {
        int nscat_step = 0;
        if (single_walk)
            nscat_step = 200;
        else
            nscat_step = 3;

        if (!scat)
            nscat_step = 0;

        for (int scat_iter = 0; scat_iter < nscat_step + 1; scat_iter++) {
            if (iso) {
                populate_spectral_flux_iso(F_down_wg_cols, // out
                                           F_up_wg_cols,   // out
                                           F_dir_wg_cols,  // in
                                           *g0_wg,         // in
                                           surface_albedo,
                                           single_walk,
                                           R_star,
                                           a,
                                           f_factor,
                                           epsi,
                                           w_0_limit,
                                           dir_beam,
                                           clouds,
                                           surface,
                                           num_cols);
            }
            else {
                populate_spectral_flux_noniso(F_down_wg_cols,
                                              F_up_wg_cols,
                                              Fc_down_wg_cols,
                                              Fc_up_wg_cols,
                                              F_dir_wg_cols,
                                              Fc_dir_wg_cols,
                                              *g0_wg_upper,
                                              *g0_wg_lower,
                                              surface_albedo,
                                              single_walk,
                                              R_star,
                                              a,
                                              f_factor,
                                              epsi,
                                              w_0_limit,
                                              delta_tau_limit,
                                              dir_beam,
                                              clouds,
                                              surface,
                                              *trans_wg_upper,
                                              *trans_wg_lower,
                                              num_cols);
            }

            cuda_check_status_or_exit(__FILE__, __LINE__);
        }

        BENCH_POINT_I_S(
            debug_nstep, debug_col_idx, "Alf_pop_spec_flx", (), ("F_up_wg", "F_down_wg"));
    }


    double* gauss_weight = *gauss_weights;
    integrate_flux(deltalambda,
                   F_down_tot_cols,
                   F_up_tot_cols,
                   F_dir_tot_cols,
                   F_net_cols,
                   F_down_wg_cols,
                   F_up_wg_cols,
                   F_dir_wg_cols,
                   F_down_band_cols,
                   F_up_band_cols,
                   F_dir_band_cols,
                   F_up_TOA_spectrum_cols,
                   gauss_weight,
                   num_cols);

    BENCH_POINT_I_S(
        debug_nstep, debug_col_idx, "Alf_int_flx", (), ("F_up_band", "F_down_band", "F_dir_band"));


    cuda_check_status_or_exit(__FILE__, __LINE__);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess) {
        printf("compute_radiative_transfer: cuda error: %s\n", cudaGetErrorString(err));
    }
}

/* Computation of temperatures, pressure and physical quantities before integrating flux */
bool alfrodull_engine::prepare_compute_flux(
    double* dev_starflux, // in: pil
    double*
                  dev_T_lay_cols, // out: it, pil, io, mmm, kil   (interpolated from T_int and then used as input to other funcs)
    double*       dev_T_int_cols,           // in: it, pii, ioi, mmmi, kii
    double*       dev_p_lay_cols,           // in: io, mmm, kil
    double*       dev_p_int_cols,           // in: ioi, mmmi, kii
    double*       dev_opac_wg_lay_cols,     // out: io
    double*       dev_opac_wg_int_cols,     // out: ioi
    double*       dev_meanmolmass_lay_cols, // out: mmm
    double*       dev_meanmolmass_int_cols, // out: mmmi
    const bool&   real_star,                // pil
    const double& fake_opac,                // io
    const bool&   interpolate_temp_and_pres,
    const bool&   interp_and_calc_flux_step,
    const int&    num_cols) {

    int nbin = opacities.nbin;

    // out: csp, cse
    int plancktable_dim  = plancktable.dim;
    int plancktable_step = plancktable.step;

    if (interpolate_temp_and_pres) {
        // it
        dim3 it_grid(int((ninterface + 15) / 16), 1, 1);
        dim3 it_block(16, 1, 1);

        interpolate_temperature<<<it_grid, it_block>>>(dev_T_lay_cols, // out
                                                       dev_T_int_cols, // in
                                                       ninterface);
        cudaDeviceSynchronize();
    }


    // pil
    dim3 pil_grid(int((nbin + 15) / 16), int(((nlayer + 2) + 15)) / 16, num_cols);
    dim3 pil_block(16, 16, 1);
    planck_interpol_layer<<<pil_grid, pil_block>>>(dev_T_lay_cols,           // in
                                                   *planckband_lay,          // out
                                                   *plancktable.planck_grid, // in
                                                   dev_starflux,             // in
                                                   null_planck_function,
                                                   real_star,
                                                   nlayer,
                                                   nbin,
                                                   plancktable_dim,
                                                   plancktable_step);
    cudaDeviceSynchronize();

    if (!iso) {
        // pii
        dim3 pii_grid(int((nbin + 15) / 16), int((ninterface + 15) / 16), num_cols);
        dim3 pii_block(16, 16, 1);
        planck_interpol_interface<<<pii_grid, pii_block>>>(dev_T_int_cols,           // in
                                                           *planckband_int,          // out
                                                           *plancktable.planck_grid, // in
                                                           null_planck_function,
                                                           ninterface,
                                                           nbin,
                                                           plancktable_dim,
                                                           plancktable_step);
        cudaDeviceSynchronize();
    }

    if (interp_and_calc_flux_step) {
        // io
        dim3 io_grid(int((nbin + 15) / 16), int((nlayer + 15) / 16), num_cols);
        dim3 io_block(16, 16, 1);
        // out -> opacities (dev_opac_wg_lay)
        // out -> scetter cross section (scatter_cross_section_...)
        interpolate_opacities<<<io_grid, io_block>>>(dev_T_lay_cols,                     // in
                                                     *opacities.dev_temperatures,        // in
                                                     dev_p_lay_cols,                     // in
                                                     *opacities.dev_pressures,           // in
                                                     *opacities.dev_kpoints,             // in
                                                     dev_opac_wg_lay_cols,               // out
                                                     *opacities.dev_scat_cross_sections, // in
                                                     *scatter_cross_section_lay,         // out
                                                     opacities.n_pressures,
                                                     opacities.n_temperatures,
                                                     opacities.ny,
                                                     nbin,
                                                     fake_opac,
                                                     nlayer + 1,
                                                     nlayer,
                                                     nlayer);


        cudaDeviceSynchronize();

        if (!iso) {

            // ioi
            dim3 ioi_grid(int((nbin + 15) / 16), int((ninterface + 15) / 16), num_cols);
            dim3 ioi_block(16, 16, 1);

            interpolate_opacities<<<ioi_grid, ioi_block>>>(dev_T_int_cols,              // in
                                                           *opacities.dev_temperatures, // in
                                                           dev_p_int_cols,              // in
                                                           *opacities.dev_pressures,    // in
                                                           *opacities.dev_kpoints,      // in
                                                           dev_opac_wg_int_cols,        // out
                                                           *opacities.dev_scat_cross_sections, // in
                                                           *scatter_cross_section_inter, // out
                                                           opacities.n_pressures,
                                                           opacities.n_temperatures,
                                                           opacities.ny,
                                                           nbin,
                                                           fake_opac,
                                                           ninterface,
                                                           ninterface,
                                                           ninterface);

            cudaDeviceSynchronize();
        }
        // mmm
        dim3 mmm_grid(int((nlayer + 15) / 16), 1, num_cols);
        dim3 mmm_block(16, 1, 1);

        meanmolmass_interpol<<<mmm_grid, mmm_block>>>(dev_T_lay_cols,              // in
                                                      *opacities.dev_temperatures, // in
                                                      dev_meanmolmass_lay_cols,    // out
                                                      *opacities.dev_meanmolmass,  // in
                                                      dev_p_lay_cols,              // in
                                                      *opacities.dev_pressures,    // in
                                                      opacities.n_pressures,
                                                      opacities.n_temperatures,
                                                      nlayer + 1,
                                                      nlayer,
                                                      nlayer);


        cudaDeviceSynchronize();

        if (!iso) {
            // mmmi
            dim3 mmmi_grid(int((ninterface + 15) / 16), 1, num_cols);
            dim3 mmmi_block(16, 1, 1);

            meanmolmass_interpol<<<mmmi_grid, mmmi_block>>>(dev_T_int_cols,              // in
                                                            *opacities.dev_temperatures, // in
                                                            dev_meanmolmass_int_cols,    // out
                                                            *opacities.dev_meanmolmass,  // in
                                                            dev_p_int_cols,              // in
                                                            *opacities.dev_pressures,    // in
                                                            opacities.n_pressures,
                                                            opacities.n_temperatures,
                                                            ninterface,
                                                            ninterface,
                                                            ninterface);


            cudaDeviceSynchronize();
        }
    }

    return true;
}

/* Integrate flux over weights and band */
void alfrodull_engine::integrate_flux(double* deltalambda,
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
                                      int     num_cols) {

    int nbin = opacities.nbin;
    int ny   = opacities.ny;

    {
        int  num_levels_per_block = 16;
        int  num_bins_per_block   = 16;
        dim3 gridsize(
            ninterface / num_levels_per_block + 1, nbin / num_bins_per_block + 1, num_cols);
        dim3 blocksize(num_levels_per_block, num_bins_per_block, 1);

        integrate_flux_band<<<gridsize, blocksize>>>(F_down_wg,
                                                     F_up_wg,
                                                     F_dir_wg,
                                                     F_down_band,
                                                     F_up_band,
                                                     F_dir_band,
                                                     F_up_TOA_spectrum,
                                                     gauss_weight,
                                                     nbin,
                                                     ninterface,
                                                     ny);

        cudaDeviceSynchronize();
    }

    {
        int  num_levels_per_block = 256;
        dim3 gridsize(ninterface / num_levels_per_block + 1, 1, num_cols);
        dim3 blocksize(num_levels_per_block, 1, 1);
        integrate_flux_tot<<<gridsize, blocksize>>>(deltalambda,
                                                    F_down_tot,
                                                    F_up_tot,
                                                    F_dir_tot,
                                                    F_net,
                                                    F_down_band,
                                                    F_up_band,
                                                    F_dir_band,
                                                    nbin,
                                                    ninterface);
        cudaDeviceSynchronize();
    }
}

/* calculate transmission function and quantities for flux computation */
void alfrodull_engine::calculate_transmission_iso(double* trans_wg,             // out
                                                  double* delta_colmass,        // in
                                                  double* opac_wg_lay,          // in
                                                  double* cloud_abs_cross_lay_, // in
                                                  double* meanmolmass_lay,      // in
                                                  double* cloud_scat_cross_lay, // in
                                                  double* g_0_cloud_lay_,       // in
                                                  double  g_0,
                                                  double  epsi,
                                                  double  epsilon2_,
                                                  double* zenith_angle_cols,
                                                  bool    scat,
                                                  bool    clouds,
                                                  int     num_cols) {

    bool hit_G_pm_denom_limit_h = false;
    // set columns wiggle iteration pingpong buffers.
    // get pointers
    // array that tells kernel to run iteration
    unsigned int* mu_star_iterate = *mu_star_iteration_buffer1;
    // array used by kernel to ask for one more iteration
    unsigned int* mu_star_iteration_request = *mu_star_iteration_buffer2;

    cudaMemset(mu_star_iterate, 1, num_cols * sizeof(unsigned int));


    int nbin = opacities.nbin;

    int ny                = opacities.ny;
    int iteration_counter = 0;


    do {
        hit_G_pm_denom_limit_h = false;
        // set global wiggle checker to 0;
        cudaMemcpy(
            *hit_G_pm_denom_limit, &hit_G_pm_denom_limit_h, sizeof(bool), cudaMemcpyHostToDevice);

        // zero out iteration request array.
        cudaMemset(mu_star_iteration_request, 0, num_cols * sizeof(unsigned int));

        dim3 grid((nbin + 15) / 16, (num_cols * ny + 3) / 4, (nlayer + 3) / 4);
        dim3 block(16, 4, 4);
        trans_iso<<<grid, block>>>(trans_wg,
                                   *delta_tau_wg,
                                   *M_term,
                                   *N_term,
                                   *P_term,
                                   *G_plus,
                                   *G_minus,
                                   delta_colmass,
                                   opac_wg_lay,
                                   cloud_abs_cross_lay_,
                                   meanmolmass_lay,
                                   *scatter_cross_section_lay,
                                   cloud_scat_cross_lay,
                                   *w0_wg,
                                   *g0_wg,
                                   g_0_cloud_lay_,
                                   g_0,
                                   epsi,
                                   epsilon2_,
                                   zenith_angle_cols,
                                   *mu_star_cols,
                                   w_0_limit,
                                   scat,
                                   nbin,
                                   ny,
                                   nlayer,
                                   num_cols,
                                   fcloud,
                                   clouds,
                                   scat_corr,
                                   mu_star_wiggle_increment,
                                   G_pm_limiter,
                                   G_pm_denom_limit,
                                   G_pm_limit_on_full_G_pm,
                                   *hit_G_pm_denom_limit,
                                   mu_star_iterate,
                                   mu_star_iteration_request,
                                   iteration_counter,
                                   debug,
                                   i2s_transition);

        cudaDeviceSynchronize();
        cudaMemcpy(
            &hit_G_pm_denom_limit_h, *hit_G_pm_denom_limit, sizeof(bool), cudaMemcpyDeviceToHost);

        // swap iteration buffer pointers
        {
            unsigned int* tmp         = mu_star_iterate;
            mu_star_iterate           = mu_star_iteration_request;
            mu_star_iteration_request = tmp;
        }
        iteration_counter += 1;

        if (hit_G_pm_denom_limit_h) {
            if (iteration_counter == wiggle_iteration_max) {
                printf("Hit maximum iteration of mu_star wiggle, bailing out\n");
                break;
            }
            printf("Hit G_pm denom limit, wiggling mu_star, loop counter %d\n", iteration_counter);
        }
    } while (hit_G_pm_denom_limit_h);
}

/* calculate transmission function and quantities for flux computation */
void alfrodull_engine::calculate_transmission_noniso(double* trans_wg_upper,
                                                     double* trans_wg_lower,
                                                     double* delta_col_upper,
                                                     double* delta_col_lower,
                                                     double* opac_wg_lay,
                                                     double* opac_wg_int,
                                                     double* cloud_abs_cross_lay_,
                                                     double* cloud_abs_cross_int_,
                                                     double* meanmolmass_lay,
                                                     double* meanmolmass_int,
                                                     double* cloud_scat_cross_lay,
                                                     double* cloud_scat_cross_int,
                                                     double* g_0_cloud_lay_,
                                                     double* g_0_cloud_int_,
                                                     double  g_0,
                                                     double  epsi,
                                                     double  epsilon2_,
                                                     double* zenith_angle_cols,
                                                     bool    scat,
                                                     bool    clouds,
                                                     int     num_cols,
                                                     int     column_idx) {


    bool hit_G_pm_denom_limit_h = false;
    // set columns wiggle iteration pingpong buffers.
    // get pointers
    // array that tells kernel to run iteration
    unsigned int* mu_star_iterate = *mu_star_iteration_buffer1;
    // array used by kernel to ask for one more iteration
    unsigned int* mu_star_iteration_request = *mu_star_iteration_buffer2;

    cudaMemset(mu_star_iterate, 1, num_cols * sizeof(unsigned int));

    int nbin = opacities.nbin;

    int ny                = opacities.ny;
    int iteration_counter = 0;

    do {
        hit_G_pm_denom_limit_h = false;
        // set global wiggle checker to 0;
        cudaMemcpy(
            *hit_G_pm_denom_limit, &hit_G_pm_denom_limit_h, sizeof(bool), cudaMemcpyHostToDevice);

        // zero out iteration request array.
        cudaMemset(mu_star_iteration_request, 0, num_cols * sizeof(unsigned int));


        dim3 grid((nbin + 15) / 16, (num_cols * ny + 3) / 4, (nlayer + 3) / 4);
        dim3 block(16, 4, 4);

        trans_noniso<<<grid, block>>>(trans_wg_upper,
                                      trans_wg_lower,
                                      *delta_tau_wg_upper,
                                      *delta_tau_wg_lower,
                                      *M_upper,
                                      *M_lower,
                                      *N_upper,
                                      *N_lower,
                                      *P_upper,
                                      *P_lower,
                                      *G_plus_upper,
                                      *G_plus_lower,
                                      *G_minus_upper,
                                      *G_minus_lower,
                                      delta_col_upper,
                                      delta_col_lower,
                                      opac_wg_lay,
                                      opac_wg_int,
                                      cloud_abs_cross_lay_,
                                      cloud_abs_cross_int_,
                                      meanmolmass_lay,
                                      meanmolmass_int,
                                      *scatter_cross_section_lay,
                                      *scatter_cross_section_inter,
                                      cloud_scat_cross_lay,
                                      cloud_scat_cross_int,
                                      *w0_wg_upper,
                                      *w0_wg_lower,
                                      *g0_wg_upper,
                                      *g0_wg_lower,
                                      g_0_cloud_lay_,
                                      g_0_cloud_int_,
                                      g_0,
                                      epsi,
                                      epsilon2_,
                                      zenith_angle_cols,
                                      *mu_star_cols,
                                      w_0_limit,
                                      scat,
                                      nbin,
                                      ny,
                                      nlayer,
                                      num_cols,
                                      fcloud,
                                      clouds,
                                      scat_corr,
                                      mu_star_wiggle_increment,
                                      G_pm_limiter,
                                      G_pm_denom_limit,
                                      G_pm_limit_on_full_G_pm,
                                      *hit_G_pm_denom_limit,
                                      mu_star_iterate,
                                      mu_star_iteration_request,
                                      iteration_counter,
                                      debug,
                                      i2s_transition,
                                      column_idx);
        cudaDeviceSynchronize();
        cudaMemcpy(
            &hit_G_pm_denom_limit_h, *hit_G_pm_denom_limit, sizeof(bool), cudaMemcpyDeviceToHost);

        // swap iteration buffer pointers
        {
            unsigned int* tmp         = mu_star_iterate;
            mu_star_iterate           = mu_star_iteration_request;
            mu_star_iteration_request = tmp;
        }
        iteration_counter += 1;

        if (hit_G_pm_denom_limit_h) {
            if (iteration_counter == wiggle_iteration_max) {
                printf("Hit maximum iteration of mu_star wiggle, bailing out\n");
                break;
            }
            printf("Hit G_pm denom limit, wiggling mu_star, loop counter %d\n", iteration_counter);
        }
    } while (hit_G_pm_denom_limit_h);
}

/* Compute direct beam from TOA spectrum, mu-star, geometry */
bool alfrodull_engine::direct_beam_flux(double* F_dir_wg,
                                        double* Fc_dir_wg,
                                        double* z_lay,
                                        double  R_planet,
                                        double  R_star,
                                        double  a,
                                        bool    dir_beam,
                                        bool    geom_zenith_corr,
                                        int     num_cols) {

    int nbin = opacities.nbin;

    int ny = opacities.ny;

    if (iso) {
        dim3 grid((ninterface + 3) / 4, (nbin + 31) / 32, (num_cols * ny + 3) / 4);
        dim3 block(4, 32, 4);
        fdir_iso<<<grid, block>>>(F_dir_wg,
                                  *planckband_lay,
                                  *delta_tau_wg,
                                  z_lay,
                                  *mu_star_cols,
                                  mu_star_limit,
                                  R_planet,
                                  R_star,
                                  a,
                                  dir_beam,
                                  geom_zenith_corr,
                                  ninterface,
                                  nbin,
                                  ny,
                                  num_cols);

        cudaDeviceSynchronize();
    }
    else {
        dim3 grid((ninterface + 3) / 4, (nbin + 31) / 32, (num_cols * ny + 3) / 4);
        dim3 block(4, 32, 4);
        fdir_noniso<<<grid, block>>>(F_dir_wg,
                                     Fc_dir_wg,
                                     *planckband_lay,
                                     *delta_tau_wg_upper,
                                     *delta_tau_wg_lower,
                                     z_lay,
                                     *mu_star_cols,
                                     mu_star_limit,
                                     R_planet,
                                     R_star,
                                     a,
                                     dir_beam,
                                     geom_zenith_corr,
                                     ninterface,
                                     nbin,
                                     ny,
                                     num_cols);

        cudaDeviceSynchronize();
    }

    return true;
}

/* Solve flux transmission matrix */
bool alfrodull_engine::populate_spectral_flux_iso_thomas(double* F_down_wg, // out
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
                                                         int     num_cols) {

    int nbin = opacities.nbin;
    int ny   = opacities.ny;


    dim3 grid((nbin + 15) / 16, (ny + 15) / 16, num_cols);
    dim3 block(16, 16, 1);
    fband_iso_thomas<<<grid, block>>>(F_down_wg,
                                      F_up_wg,
                                      F_dir_wg,
                                      *planckband_lay,
                                      *w0_wg,
                                      *M_term,
                                      *N_term,
                                      *P_term,
                                      *G_plus,
                                      *G_minus,
                                      *A_buff,       // thomas worker
                                      *B_buff,       // thomas worker
                                      *C_buff,       // thomas worker
                                      *D_buff,       // thomas worker
                                      *C_prime_buff, // thomas worker
                                      *D_prime_buff, // thomas worker
                                      *X_buff,       // thomas worker
                                      g_0_tot,
                                      surface_albedo,
                                      singlewalk,
                                      Rstar,
                                      a,
                                      ninterface,
                                      nbin,
                                      f_factor,
                                      *mu_star_cols,
                                      ny,
                                      num_cols,
                                      epsi,
                                      dir_beam,
                                      clouds,
                                      scat_corr,
                                      debug,
                                      i2s_transition);

    cudaDeviceSynchronize();

    return true;
}

bool alfrodull_engine::populate_spectral_flux_iso(double* F_down_wg, // out
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
                                                  int     num_cols) {

    int nbin = opacities.nbin;
    int ny   = opacities.ny;

    dim3 grid((nbin + 15) / 16, (ny + 15) / 16, num_cols);
    dim3 block(16, 16, 1);
    fband_iso_notabu<<<grid, block>>>(F_down_wg,
                                      F_up_wg,
                                      F_dir_wg,
                                      *planckband_lay,
                                      *w0_wg,
                                      *M_term,
                                      *N_term,
                                      *P_term,
                                      *G_plus,
                                      *G_minus,
                                      g_0_tot,
                                      surface_albedo,
                                      singlewalk,
                                      Rstar,
                                      a,
                                      ninterface,
                                      nbin,
                                      f_factor,
                                      *mu_star_cols,
                                      ny,
                                      num_cols,
                                      epsi,
                                      dir_beam,
                                      clouds,
                                      scat_corr,
                                      debug,
                                      i2s_transition);

    cudaDeviceSynchronize();
    return true;
}

// calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
bool alfrodull_engine::populate_spectral_flux_noniso(double* F_down_wg,
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
                                                     int     num_cols) {
    int nbin = opacities.nbin;
    int ny   = opacities.ny;


    dim3 grid((nbin + 15) / 16, (ny + 15) / 16, num_cols);
    dim3 block(16, 16, 1);
    // calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
    fband_noniso_notabu<<<grid, block>>>(F_down_wg,
                                         F_up_wg,
                                         Fc_down_wg,
                                         Fc_up_wg,
                                         F_dir_wg,
                                         Fc_dir_wg,
                                         *planckband_lay,
                                         *planckband_int,
                                         *w0_wg_upper,
                                         *w0_wg_lower,
                                         *delta_tau_wg_upper,
                                         *delta_tau_wg_lower,
                                         *M_upper,
                                         *M_lower,
                                         *N_upper,
                                         *N_lower,
                                         *P_upper,
                                         *P_lower,
                                         *G_plus_upper,
                                         *G_plus_lower,
                                         *G_minus_upper,
                                         *G_minus_lower,
                                         g_0_tot_upper,
                                         g_0_tot_lower,
                                         surface_albedo,
                                         singlewalk,
                                         Rstar,
                                         a,
                                         ninterface,
                                         nbin,
                                         f_factor,
                                         *mu_star_cols,
                                         ny,
                                         num_cols,
                                         epsi,
                                         delta_tau_limit,
                                         dir_beam,
                                         clouds,
                                         scat_corr,
                                         debug,
                                         i2s_transition);

    cudaDeviceSynchronize();

    return true;
}

// calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
bool alfrodull_engine::populate_spectral_flux_noniso_thomas(double* F_down_wg,
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
                                                            int     num_cols) {
    int nbin = opacities.nbin;
    int ny   = opacities.ny;

    dim3 block(16, 16, 1);

    dim3 grid((nbin + 15) / 16, (ny + 15) / 16, num_cols);

    // calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
    fband_noniso_thomas<<<grid, block>>>(F_down_wg,
                                         F_up_wg,
                                         F_dir_wg,
                                         Fc_dir_wg,
                                         *planckband_lay,
                                         *planckband_int,
                                         *w0_wg_upper,
                                         *w0_wg_lower,
                                         *delta_tau_wg_upper,
                                         *delta_tau_wg_lower,
                                         *M_upper,
                                         *M_lower,
                                         *N_upper,
                                         *N_lower,
                                         *P_upper,
                                         *P_lower,
                                         *G_plus_upper,
                                         *G_plus_lower,
                                         *G_minus_upper,
                                         *G_minus_lower,
                                         *A_buff,       // thomas worker
                                         *B_buff,       // thomas worker
                                         *C_buff,       // thomas worker
                                         *D_buff,       // thomas worker
                                         *C_prime_buff, // thomas worker
                                         *D_prime_buff, // thomas worker
                                         *X_buff,       // thomas worker
                                         g_0_tot_upper,
                                         g_0_tot_lower,
                                         surface_albedo,
                                         singlewalk,
                                         Rstar,
                                         a,
                                         ninterface,
                                         nbin,
                                         f_factor,
                                         *mu_star_cols,
                                         ny,
                                         1,
                                         epsi,
                                         delta_tau_limit,
                                         dir_beam,
                                         clouds,
                                         scat_corr,
                                         debug,
                                         i2s_transition);

    cudaDeviceSynchronize();

    return true;
}

/* Compute g0 and w0 values for output */
bool alfrodull_engine::get_column_integrated_g0_w0(double* g0_, double* w0_, const int& num_cols) {

    int nbin = opacities.nbin;
    int ny   = opacities.ny;

    if (!iso) {
        // compute mean of upper and lower band
        int  num_val              = nlayer * nbin * ny * num_cols;
        int  num_levels_per_block = 256;
        dim3 gridsize(num_val / num_levels_per_block + 1);
        dim3 blocksize(num_levels_per_block);

        arrays_mean<<<gridsize, blocksize>>>(*w0_wg_upper, *w0_wg_lower, *w0_wg, num_val);
        arrays_mean<<<gridsize, blocksize>>>(*g0_wg_upper, *g0_wg_lower, *g0_wg, num_val);
    }

    {
        int  num_levels_per_block = 16;
        int  num_bins_per_block   = 16;
        dim3 gridsize(nlayer / num_levels_per_block + 1, nbin / num_bins_per_block + 1, num_cols);
        dim3 blocksize(num_levels_per_block, num_bins_per_block, 1);

        integrate_val_band<<<gridsize, blocksize>>>(*w0_wg, w0_, *gauss_weights, nbin, nlayer, ny);
        integrate_val_band<<<gridsize, blocksize>>>(*g0_wg, g0_, *gauss_weights, nbin, nlayer, ny);

        cudaDeviceSynchronize();
    }

    // This would integrate over bands
    // {
    //     int  num_levels_per_block = 256;
    //     dim3 gridsize(ninterface / num_levels_per_block + 1);
    //     dim3 blocksize(num_levels_per_block);

    //     double* deltalambda = *opacities.dev_opac_deltawave;

    //     integrate_val_tot<<<gridsize, blocksize>>>(g0_, *g0_band, deltalambda, nbin, nlayer);
    //     integrate_val_tot<<<gridsize, blocksize>>>(w0_, *w0_band, deltalambda, nbin, nlayer);

    //     cudaDeviceSynchronize();
    // }
    return true;
}
