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
// Two stream radiative transfer class to use in THOR
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

#include "two_streams_radiative_transfer.h"

#include "binary_test.h"
#include "debug.h"
#include "debug_helpers.h"

#include "alfrodull_engine.h"

#include "physics_constants.h"

#include "directories.h"
#include "storage.h"

#include <string>

#include <functional>
#include <map>

#include "insolation.h"

#include "math_helpers.h"

USE_BENCHMARK();

using std::string;

// show progress bar
// not really necessary when running on multiple columns at low res
//#define COLUMN_LOOP_PROGRESS_BAR

// debugging printout
//#define DEBUG_PRINTOUT_ARRAYS
// dump TP profile to run in HELIOS for profile comparison
//#define DUMP_HELIOS_TP
// stride for column TP profile dump
#ifdef DUMP_HELIOS_TP
const int HELIOS_TP_STRIDE = 1;
#endif // DUMP_HELIOS_TP

//***************************************************************************************************

const char PBSTR[] = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
const int  PBWIDTH = 60;

void print_progress(double percentage) {
    int val  = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

two_streams_radiative_transfer::two_streams_radiative_transfer() {
}

two_streams_radiative_transfer::~two_streams_radiative_transfer() {
}

void two_streams_radiative_transfer::print_config() {
    log::printf("    Stellar Temperature:                     %g K\n", T_star);
    log::printf("    Internal Temperature:                    %g K\n", T_internal);
    log::printf("    Stellar Radius:                          %g [R_SUN]\n", R_star_config);
    log::printf("    Planet-Star Distance:                    %g [au]\n", planet_star_dist_config);

    log::printf("    Real Stellar Spectrum:                   %s\n", real_star ? "true" : "false");
    log::printf("    Stellar Spectrum File:                   %s\n", stellar_spectrum_file.c_str());

    log::printf("    Surface Albedo File:                     %s\n", surface_albedo_file.c_str());

    log::printf("    Isotermal Single Layers:                 %s\n", iso ? "true" : "false");

    log::printf("    Thomas Solver:                           %s\n", thomas ? "true" : "false");
    log::printf("    Iterative Solver use high def:           %s\n",
                scat_single_walk ? "true" : "false");
    log::printf("    Experimental Constant Opacity Offsets:   %g\n", experimental_opacities_offset);

    log::printf("    Gas constant g0 (without clouds):        %g\n", g_0);
    log::printf("    epsilon:                                 %g\n", epsi);
    log::printf("    epsilon_2:                               %g\n", epsilon_2);
    log::printf("    Opacity Cutoff:                          %g\n", fake_opac);

    log::printf("    Opacities File:                          %s\n", opacities_file.c_str());
    log::printf("    Opacities File is CGS:                   %s\n",
                opacity_file_is_CGS ? "true" : "false");

    log::printf("    Scattering:                              %s\n", scat ? "true" : "false");
    log::printf("    Apply Scattering Correction (else E=1):  %s\n", scat_corr ? "true" : "false");

    log::printf("    Clouds:                                  %s\n", clouds ? "true" : "false");
    log::printf("    fcloud:                                  %g\n", fcloud);
    log::printf("    Cloud Properties File:                   %s\n", cloud_filename.c_str());
    log::printf("    Store w0 g0 (per band):                  %s\n",
                store_w0_g0 ? "true" : "false");

    log::printf("    Store dir flux spectrum (per band):      %s\n",
                store_dir_spectrum ? "true" : "false");

    log::printf("    Null Planck Function:                    %s\n",
                null_planck_function ? "true" : "false");

    log::printf("    Direct Beam:                             %s\n", dir_beam ? "true" : "false");
    // log::printf("    Geometrical Zenith Correction:           %s\n",
    //             geom_zenith_corr ? "true" : "false");
    log::printf("    Direct Beam Tangent Angle Limit:         %g°\n", mu_star_limit_degrees);

    log::printf("    w0 limit:                                %g\n", w_0_limit);
    log::printf("    i2s transition:                          %g\n", i2s_transition);

    log::printf("    Compute Every N step:                    %d\n", compute_every_n_iteration);
    log::printf("    Number of Parallel Columns:              %d\n", num_parallel_columns);

    // log::printf("    Apply G_pm limiter:                      %s\n",
    //             G_pm_limiter ? "true" : "false");
    // log::printf("    G_pm limit on full G_pm:                 %g\n", G_pm_limit_on_full_G_pm);
    // log::printf("    G_pm limit:                              %g\n", G_pm_denom_limit);
    // log::printf("    G_pm angle increment:                    %g\n", mu_star_wiggle_increment);
    // log::printf("    mu_star wiggle max iterations:           %d\n", wiggle_iteration_max);

    // spinup-spindown parameters
    log::printf("    Spin up start step:                      %d\n", spinup_start_step);
    log::printf("    Spin up stop step:                       %d\n", spinup_stop_step);
    log::printf("    Spin down start step:                    %d\n", spindown_start_step);
    log::printf("    Spin down stop step:                     %d\n", spindown_stop_step);
    log::printf("    Debug output:                            %s\n",
                debug_output ? "true" : "false");
}

bool two_streams_radiative_transfer::configure(config_file &config_reader) {
    // variables reused from DG
    config_reader.append_config_var("Tstar", T_star, T_star);
    config_reader.append_config_var("Tint", T_internal, T_internal);
    config_reader.append_config_var(
        "planet_star_dist", planet_star_dist_config, planet_star_dist_config);
    config_reader.append_config_var("radius_star", R_star_config, R_star_config);

    config_reader.append_config_var("Alf_thomas", thomas, thomas);
    config_reader.append_config_var("Alf_scat_single_walk", scat_single_walk, scat_single_walk);
    config_reader.append_config_var(
        "Alf_exp_opac_offset", experimental_opacities_offset, experimental_opacities_offset);
    config_reader.append_config_var("Alf_iso", iso, iso);
    config_reader.append_config_var("Alf_real_star", real_star, real_star);
    config_reader.append_config_var(
        "Alf_stellar_spectrum", stellar_spectrum_file, stellar_spectrum_file);
    config_reader.append_config_var("Alf_surface_albedo", surface_albedo_file, surface_albedo_file);
    config_reader.append_config_var("Alf_fake_opac", fake_opac, fake_opac);

    config_reader.append_config_var("Alf_g_0", g_0, g_0);
    config_reader.append_config_var("Alf_espilon_2", epsilon_2, epsilon_2);
    // config_reader.append_config_var("Alf_G_pm_max_limiter", G_pm_limiter, G_pm_limiter);
    // config_reader.append_config_var("Alf_G_pm_limit", G_pm_denom_limit, G_pm_denom_limit);
    // config_reader.append_config_var(
    //     "Alf_G_pm_limit_on_G_pm", G_pm_limit_on_full_G_pm, G_pm_limit_on_full_G_pm);
    // config_reader.append_config_var(
    //     "Alf_G_pm_mu_star_increment", mu_star_wiggle_increment, mu_star_wiggle_increment);
    // config_reader.append_config_var(
    //     "Alf_mu_star_iteration_max", wiggle_iteration_max, wiggle_iteration_max);
    config_reader.append_config_var(
        "Alf_direct_beam_angle_limit", mu_star_limit_degrees, mu_star_limit_degrees);
    config_reader.append_config_var("Alf_scat", scat, scat);
    config_reader.append_config_var("Alf_scat_corr", scat_corr, scat_corr);

    config_reader.append_config_var("Alf_dir_beam", dir_beam, dir_beam);
    //config_reader.append_config_var("Alf_geom_zenith_corr", geom_zenith_corr, geom_zenith_corr);
    config_reader.append_config_var("Alf_i2s_transition", i2s_transition, i2s_transition);

    config_reader.append_config_var("Alf_opacities_file", opacities_file, opacities_file);
    config_reader.append_config_var(
        "Alf_opacities_file_in_CGS", opacity_file_is_CGS, opacity_file_is_CGS);
    config_reader.append_config_var(
        "Alf_compute_every_nstep", compute_every_n_iteration, compute_every_n_iteration);

    // spin up spin down
    config_reader.append_config_var("Alf_spinup_start", spinup_start_step, spinup_start_step);
    config_reader.append_config_var("Alf_spinup_stop", spinup_stop_step, spinup_stop_step);
    config_reader.append_config_var("Alf_spindown_start", spindown_start_step, spindown_start_step);
    config_reader.append_config_var("Alf_spindown_stop", spindown_stop_step, spindown_stop_step);

    config_reader.append_config_var(
        "Alf_num_parallel_columns", num_parallel_columns, num_parallel_columns);

    config_reader.append_config_var("Alf_clouds", clouds, clouds);
    config_reader.append_config_var("Alf_fcloud", fcloud, fcloud);
    config_reader.append_config_var("Alf_cloudfile", cloud_filename, cloud_filename);

    config_reader.append_config_var("Alf_store_w0_g0", store_w0_g0, store_w0_g0);
    config_reader.append_config_var(
        "Alf_store_dir_spectrum", store_dir_spectrum, store_dir_spectrum);
    config_reader.append_config_var(
        "Alf_null_planck_function", null_planck_function, null_planck_function);


    config_reader.append_config_var("Alf_debug", debug_output, debug_output);

    return true;
}

bool two_streams_radiative_transfer::initialise_memory(
    const ESP &              esp,
    device_RK_array_manager &phy_modules_core_arrays) {
    bool out = true;
    nlayer   = esp.nv;

    R_star_SI = R_star_config * R_SUN;

    planet_star_dist_SI = planet_star_dist_config * AU;

    if (num_parallel_columns < 1)
        num_parallel_columns = 1;
    // as set in host_functions.set_up_numerical_parameters
    // w_0_limit
    w_0_limit = 1.0 - 1e-14;

    double f_factor = 1.0;

    epsi = 1.0 / diffusivity;

    alf.thomas = thomas;

    alf.G_pm_limiter             = G_pm_limiter;
    alf.G_pm_denom_limit         = G_pm_denom_limit;
    alf.mu_star_wiggle_increment = mu_star_wiggle_increment;

    double mu_star_limit = cos((90.0 + mu_star_limit_degrees) / 180.0 * M_PI);

    alf.set_parameters(nlayer,    // const int&    nlayer_,
                       iso,       // const bool&   iso_,
                       T_star,    // const double& T_star_,
                       real_star, // const bool&   real_star_,
                       null_planck_function,
                       fake_opac,           // const double& fake_opac_,
                       g_0,                 // const double& g_0_,
                       epsi,                // const double& epsi_,
                       epsilon_2,           // const double& epsilon_2_,
                       scat,                // const bool&   scat_,
                       scat_corr,           // const bool&   scat_corr_,
                       0.0,                 // const double& R_planet_, filled in later
                       R_star_SI,           // const double& R_star_,
                       planet_star_dist_SI, // const double& a_,
                       dir_beam,            // const bool&   dir_beam_,
                       //geom_zenith_corr,    // const bool&   geom_zenith_corr_,
                       f_factor,       // const double& f_factor_,
                       w_0_limit,      // const double& w_0_limit_,
                       i2s_transition, // const double& i2s_transition_,
                       mu_star_limit,
                       wiggle_iteration_max,
                       G_pm_limit_on_full_G_pm,
                       num_parallel_columns,
                       debug_output); // const bool&   debug_

    // initialise opacities table -> gives frequency bins
    // set opacity offset for test
    alf.set_experimental_opacity_offset(experimental_opacities_offset);


    if (!path_exists(opacities_file)) {
        log::printf("Opacity file not found: %s\n", opacities_file.c_str());
        exit(EXIT_FAILURE);
    }

    alf.load_opacities(opacities_file, opacity_file_is_CGS);

    cudaDeviceSynchronize();
    log::printf("Loaded opacities, using %d bins with %d weights per bin\n",
                alf.opacities.nbin,
                alf.opacities.ny);

    alf.allocate_internal_variables();
    cuda_check_status_or_exit(__FILE__, __LINE__);
    
    int ninterface         = nlayer + 1;
    int nlayer_plus1       = nlayer + 1;
    int nbin               = alf.opacities.nbin;
    int ny                 = alf.opacities.ny;
    int nlayer_nbin        = nlayer * nbin;
    int ninterface_nbin    = ninterface * nbin;
    int ninterface_wg_nbin = ninterface * ny * nbin;
    int ncol               = num_parallel_columns;

    if (real_star) {
        // load star flux.
        std::printf("Using Stellar Flux file %s\n", stellar_spectrum_file.c_str());

        if (!path_exists(stellar_spectrum_file)) {
            log::printf("Stellar spectrum file not found: %s\n", stellar_spectrum_file.c_str());
            exit(EXIT_FAILURE);
        }
        
        star_flux.allocate(nbin);
        
        double lambda_spectrum_scale = 1.0;
        double flux_scale            = 1.0;

        storage s(stellar_spectrum_file, true);
        if (s.has_table("wavelength") && s.has_table("flux")) {
            std::unique_ptr<double[]> lambda_ptr  = nullptr;
            int                       lambda_size = 0;

            std::unique_ptr<double[]> flux_ptr  = nullptr;
            int                       flux_size = 0;

            s.read_table("wavelength", lambda_ptr, lambda_size);
            s.read_table("flux", flux_ptr, flux_size);

            if (lambda_size != nbin || lambda_size != flux_size) {
                log::printf("Wrong size for stellar size arrays\n");
                log::printf("Lambda: %d\n", lambda_size);
                log::printf("Flux: %d\n", flux_size);
                log::printf("nbin: %d\n", nbin);
                exit(EXIT_FAILURE);
            }

            bool                      lambda_check = true;
            double                    epsilon      = 1e-4;
            std::shared_ptr<double[]> star_flux_h  = star_flux.get_host_data_ptr();
            for (int i = 0; i < nbin; i++) {
                star_flux_h[i] = flux_ptr[i] * flux_scale;
                bool check =
                    fabs(lambda_ptr[i] * lambda_spectrum_scale - alf.opacities.data_opac_wave[i])
                        / alf.opacities.data_opac_wave[i]
                    < epsilon;

                if (!check)
                    printf("Missmatch in wavelength at idx [%d] l_spectrum(%g) != "
                           "l_opac(%g) \n",
                           i,
                           lambda_ptr[i] * lambda_spectrum_scale,
                           alf.opacities.data_opac_wave[i]);
                lambda_check &= check;
            }

            star_flux.put();

            if (!lambda_check) {
                log::printf("wavelength points mismatch between stellar spectrum and "
                            "opacities\n");
                exit(EXIT_FAILURE);
            }
        }
        else {
            log::printf("table wavelength or flux not found in stellar flux file\n");
            exit(EXIT_FAILURE);
        }
        printf("Stellar flux loaded\n");
    }

    cuda_check_status_or_exit(__FILE__, __LINE__);

    if (esp.surface) {
        //read in surface albedo file
        // load star flux.
        std::printf("Using surface albedo file %s\n", surface_albedo_file.c_str());
        surface_albedo.allocate(nbin);
        if (!path_exists(surface_albedo_file)) {
            log::printf("Surface albedo file not found: %s\n", surface_albedo_file.c_str());
            exit(EXIT_FAILURE);
        }

        double lambda_spectrum_scale = 1.0;
        //double flux_scale            = 1.0;

        storage s(surface_albedo_file, true);
        if (s.has_table("wavelength") && s.has_table("albedo")) {
            std::unique_ptr<double[]> lambda_ptr  = nullptr;
            int                       lambda_size = 0;

            std::unique_ptr<double[]> albedo_ptr  = nullptr;
            int                       albedo_size = 0;

            s.read_table("wavelength", lambda_ptr, lambda_size);
            s.read_table("albedo", albedo_ptr, albedo_size);

            if (lambda_size != nbin || lambda_size != albedo_size) {
                log::printf("Wrong size for albedo size arrays\n");
                log::printf("Lambda: %d\n", lambda_size);
                log::printf("Albedo: %d\n", albedo_size);
                log::printf("nbin: %d\n", nbin);
                exit(EXIT_FAILURE);
            }

            bool                      lambda_check     = true;
            double                    epsilon          = 1e-4;
            std::shared_ptr<double[]> surface_albedo_h = surface_albedo.get_host_data_ptr();
            for (int i = 0; i < nbin; i++) {
                surface_albedo_h[i] = albedo_ptr[i];
                bool check =
                    fabs(lambda_ptr[i] * lambda_spectrum_scale - alf.opacities.data_opac_wave[i])
                        / alf.opacities.data_opac_wave[i]
                    < epsilon;

                if (!check)
                    printf("Missmatch in wavelength at idx [%d] l_albedo(%g) != "
                           "l_opac(%g) \n",
                           i,
                           lambda_ptr[i] * lambda_spectrum_scale,
                           alf.opacities.data_opac_wave[i]);
                lambda_check &= check;
            }

            surface_albedo.put();

            if (!lambda_check) {
                log::printf("wavelength points mismatch between surface albedo and "
                            "opacities\n");
                exit(EXIT_FAILURE);
            }
        }
        else {
            log::printf("table wavelength or flux not found in surface albedo file\n");
            exit(EXIT_FAILURE);
        }
        printf("Surface albedo loaded\n");
    }

    else { // surface off, set albedo to 1.0 everywhere
        surface_albedo.allocate(nbin);
        std::shared_ptr<double[]> surface_albedo_h = surface_albedo.get_host_data_ptr();
        for (int i = 0; i < nbin; i++) {
            surface_albedo_h[i] = 1.0;
        }

        surface_albedo.put();
    }
    // allocate interface state variables to be interpolated

    pressure_int.allocate(ncol * ninterface);
    temperature_int.allocate(ncol * ninterface);
    temperature_lay.allocate(ncol * nlayer_plus1);
    density_int.allocate(ncol * ninterface);

    F_down_wg.allocate(ncol * ninterface_wg_nbin);
    F_up_wg.allocate(ncol * ninterface_wg_nbin);
    F_dir_wg.allocate(ncol * ninterface_wg_nbin);

    cuda_check_status_or_exit(__FILE__, __LINE__);

    if (!iso) {
        if (!thomas) {
            Fc_down_wg.allocate(ncol * ninterface_wg_nbin);
            Fc_up_wg.allocate(ncol * ninterface_wg_nbin);
        }
        Fc_dir_wg.allocate(ncol * ninterface_wg_nbin);
    }

    F_down_tot.allocate(esp.point_num * ninterface);
    F_up_tot.allocate(esp.point_num * ninterface);
    F_dir_tot.allocate(esp.point_num * ninterface);
    F_down_band.allocate(ncol * ninterface_nbin);
    F_up_band.allocate(ncol * ninterface_nbin);
    if (store_dir_spectrum)
        F_dir_band.allocate(esp.point_num * ninterface_nbin);
    else
        F_dir_band.allocate(ncol * ninterface_nbin);
    F_net.allocate(esp.point_num * ninterface);

    F_up_TOA_spectrum.allocate(esp.point_num * nbin);

    Qheat.allocate(esp.point_num * nlayer);

    cuda_check_status_or_exit(__FILE__, __LINE__);

    if (store_w0_g0) {
        // output for storage
        g0_tot.allocate(esp.point_num * nlayer_nbin);
        w0_tot.allocate(esp.point_num * nlayer_nbin);
        // output for storage
        g0_tot.zero();
        w0_tot.zero();
    }

    cuda_check_status_or_exit(__FILE__, __LINE__);

    if (clouds) {
        // load cloud file
        if (!path_exists(cloud_filename)) {
            log::printf("Cloud file not found: %s\n", cloud_filename.c_str());
            exit(EXIT_FAILURE);
        }
      

        alf.cloud_opacities.load(cloud_filename);

        int  asymmetry_size     = alf.cloud_opacities.dev_asymmetry.get_size();
        int  scat_cs_size       = alf.cloud_opacities.dev_scat_cross_sections.get_size();
        int  abs_cs_size        = alf.cloud_opacities.dev_abs_cross_sections.get_size();
        bool cloud_load_failure = false;

        std::shared_ptr<double[]> cloud_wavelength_h =
            alf.cloud_opacities.dev_wavelength.get_host_data();

        if (asymmetry_size != nbin) {
            log::printf("Wrong size for cloud assymetry array size\n");
            log::printf("size: %d, nbin: %d\n", asymmetry_size, nbin);
            cloud_load_failure = true;
        }

        if (scat_cs_size != nbin) {
            log::printf("Wrong size for cloud scattering cross section array size\n");
            log::printf("size: %d, nbin: %d\n", scat_cs_size, nbin);
            cloud_load_failure = true;
        }

        if (abs_cs_size != nbin) {
            log::printf("Wrong size for cloud absorption cross section array size\n");
            log::printf("size: %d, nbin: %d\n", abs_cs_size, nbin);
            cloud_load_failure = true;
        }

        double lambda_spectrum_scale = 1.0;

        bool   lambda_check = true;
        double epsilon      = 1e-4;
        for (int i = 0; i < nbin; i++) {
            bool check = fabs(cloud_wavelength_h[i] * lambda_spectrum_scale
                              - alf.opacities.data_opac_wave[i])
                             / alf.opacities.data_opac_wave[i]
                         < epsilon;

            if (!check)
                printf("Missmatch in wavelength at idx [%d] l_cloud(%g) != "
                       "l_opac(%g) \n",
                       i,
                       cloud_wavelength_h[i] * lambda_spectrum_scale,
                       alf.opacities.data_opac_wave[i]);
            lambda_check &= check;
        }

        if (!lambda_check) {
            log::printf("wavelength points mismatch between cloud spectrum and "
                        "opacities\n");
            exit(EXIT_FAILURE);
        }

        if (cloud_load_failure) {
            exit(EXIT_FAILURE);
        }

        alf.set_clouds_data(clouds,
                            *alf.cloud_opacities.dev_abs_cross_sections,
                            *alf.cloud_opacities.dev_abs_cross_sections,
                            *alf.cloud_opacities.dev_scat_cross_sections,
                            *alf.cloud_opacities.dev_scat_cross_sections,
                            *alf.cloud_opacities.dev_asymmetry,
                            *alf.cloud_opacities.dev_asymmetry,
                            fcloud);
    }
    else {
        // all clouds set to zero. Not used.
        g_0_tot_lay.allocate(nbin);
        g_0_tot_int.allocate(nbin);
        cloud_abs_cross_lay.allocate(nbin);
        cloud_abs_cross_int.allocate(nbin);
        cloud_scat_cross_lay.allocate(nbin);
        cloud_scat_cross_int.allocate(nbin);

        g_0_tot_lay.zero();
        g_0_tot_int.zero();
        cloud_abs_cross_lay.zero();
        cloud_abs_cross_int.zero();
        cloud_scat_cross_lay.zero();
        cloud_scat_cross_int.zero();

        fcloud = 0.0;

        alf.set_clouds_data(clouds,
                            *cloud_abs_cross_lay,
                            *cloud_abs_cross_int,
                            *cloud_scat_cross_lay,
                            *cloud_scat_cross_int,
                            *g_0_tot_lay,
                            *g_0_tot_int,
                            fcloud);

        cuda_check_status_or_exit(__FILE__, __LINE__);
    }

    cudaError_t err = cudaGetLastError();

    // Check device query
    if (err != cudaSuccess) {
        log::printf("[%s:%d] CUDA error check reports error: %s\n",
                    __FILE__,
                    __LINE__,
                    cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

#ifdef BENCHMARKING
    std::map<string, output_def> debug_arrays = {
        {"F_net", {F_net.ptr_ref(), esp.point_num * ninterface, "Fnet", "Fn", true, dummy}},

        {"F_up_tot",
         {F_up_tot.ptr_ref(), esp.point_num * ninterface, "Fuptot", "Fut", true, dummy}},
        {"F_down_tot",
         {F_down_tot.ptr_ref(), esp.point_num * ninterface, "Fdowntot", "Fdt", true, dummy}},
        {"F_up_wg", {F_up_wg.ptr_ref(), ncol * ninterface_wg_nbin, "Fupwg", "Fuw", true, dummy}},
        {"F_down_wg",
         {F_down_wg.ptr_ref(), ncol * ninterface_wg_nbin, "Fdownwg", "Fdw", true, dummy}},
        {"F_up_band", {F_up_band.ptr_ref(), ncol * ninterface_nbin, "Fupband", "Fub", true, dummy}},
        {"F_down_band",
         {F_down_band.ptr_ref(), ncol * ninterface_nbin, "Fdownband", "Fdb", true, dummy}},
        {"F_dir_wg",
         {F_dir_wg.ptr_ref(), ncol * ninterface_wg_nbin, "Fdirwg", "Fdirw", true, dummy}},

        {"F_dir_band",
         {F_dir_band.ptr_ref(), ncol * ninterface_nbin, "Fdirband", "Fdib", true, dummy}},

        {"T_lay", {temperature_lay.ptr_ref(), ncol * nlayer_plus1, "T_lay", "Tl", true, dummy}},
        {"T_int", {temperature_int.ptr_ref(), ncol * ninterface, "T_int", "Ti", true, dummy}},
        {"P_int", {pressure_int.ptr_ref(), ncol * ninterface, "P_int", "Pi", true, dummy}},

        //        {"col_mu_star", {col_mu_star.ptr_ref(), esp.point_num,
        //        "col_mu_star", "cMu", true, dummy}},
        {"AlfQheat", {Qheat.ptr_ref(), esp.point_num * nlayer, "AlfQheat", "aQh", true, dummy}}};
    BENCH_POINT_REGISTER_PHY_VARS(debug_arrays, (), ());
#endif // BENCHMARKING
    return out;
}

bool two_streams_radiative_transfer::initial_conditions(const ESP &            esp,
                                                        const SimulationSetup &sim,
                                                        storage *              s) {
    if (spinup_start_step > -1 || spinup_stop_step > -1) {
        if (spinup_stop_step < spinup_start_step)
            printf("Alf: inconsistent spinup_start (%d) and spinup_stop (%d) values\n",
                   spinup_start_step,
                   spinup_stop_step);
    }
    if (spindown_start_step > -1 || spindown_stop_step > -1) {
        if (spindown_stop_step < spindown_start_step)
            printf("Alf: inconsistent spindown_start (%d) and spindown_stop (%d) "
                   "values\n",
                   spindown_start_step,
                   spindown_stop_step);
    }

    bool out = true;
    // what should be initialised here and what is to initialise at each loop ?
    // what to initialise here and what to do in initialise memory ?

    // this is only known here, comes from sim setup.
    alf.R_planet = sim.A;
    cuda_check_status_or_exit(__FILE__, __LINE__);
    // initialise planck tables
    alf.prepare_planck_table();
    log::printf("Built Planck Table for %d bins, Star temp %g K\n", alf.opacities.nbin, alf.T_star);
    // initialise alf

    alf.correct_incident_energy(*star_flux, real_star, true);

    // internal flux from internal temperature
    F_intern = STEFANBOLTZMANN * pow(T_internal, 4);

    cuda_check_status_or_exit(__FILE__, __LINE__);

    // request insolation computation
    esp.insolation.set_require();
    return out;
}

// ************************************************************************************************
// initialise delta_colmass arrays from pressure
// same as helios.source.host_functions.construct_grid
__global__ void initialise_delta_colmass_pressure_noniso(double *delta_col_mass_upper_cols,
                                                         double *delta_col_mass_lower_cols,
                                                         double *pressure_lay_cols,
                                                         double *pressure_int_cols,
                                                         double  gravit,
                                                         int     num_layers,
                                                         int     num_columns) {
    int layer_idx = blockIdx.x * blockDim.x + threadIdx.x;
    //    int col_idx   = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;

    if (layer_idx < num_layers) {
        // get offset into start of column data
        double *delta_col_mass_upper = &(delta_col_mass_upper_cols[col_block_idx * num_layers]);
        double *delta_col_mass_lower = &(delta_col_mass_lower_cols[col_block_idx * num_layers]);
        double *pressure_int         = &(pressure_int_cols[col_block_idx * (num_layers + 1)]);
        double *pressure_lay         = &(pressure_lay_cols[col_block_idx * num_layers]);

        delta_col_mass_upper[layer_idx] =
            fabs(pressure_lay[layer_idx] - pressure_int[layer_idx + 1]) / gravit;
        delta_col_mass_lower[layer_idx] =
            fabs(pressure_int[layer_idx] - pressure_lay[layer_idx]) / gravit;

        if ((delta_col_mass_upper[layer_idx] < 0.0) || (delta_col_mass_lower[layer_idx] < 0.0))
            printf("Negative delta_col_mass (%g, %g), layer: %d, col: %d\n",
                   delta_col_mass_upper[layer_idx],
                   delta_col_mass_lower[layer_idx],
                   layer_idx,
                   col_block_idx);
    }
}

// initialise delta_colmass arrays from pressure
// same as helios.source.host_functions.construct_grid
__global__ void initialise_delta_colmass_pressure_iso(double *delta_col_mass_cols,
                                                      double *pressure_int_cols,
                                                      double  gravit,
                                                      int     num_layers,
                                                      int     num_columns) {
    int layer_idx = blockIdx.x * blockDim.x + threadIdx.x;
    // int col_idx   = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;
    int col_size      = num_layers;

    if (layer_idx < num_layers) {
        // get offset into start of column data
        double *delta_col_mass = &(delta_col_mass_cols[col_block_idx * col_size]);
        double *pressure_int   = &(pressure_int_cols[col_block_idx * (num_layers + 1)]);
        delta_col_mass[layer_idx] =
            (pressure_int[layer_idx] - pressure_int[layer_idx + 1]) / gravit;
        if (delta_col_mass[layer_idx] < 0.0)
            printf("Negative delta_col_mass (%g), layer: %d, col: %d\n",
                   delta_col_mass[layer_idx],
                   layer_idx,
                   col_block_idx);
    }
}

// initialise delta_colmass arrays from pressure
// same as helios.source.host_functions.construct_grid
__global__ void initialise_delta_colmass_density_noniso(double *delta_col_mass_upper_cols,
                                                        double *delta_col_mass_lower_cols,
                                                        double *altitude_lay,
                                                        double *altitude_int,
                                                        double *density_lay_cols,
                                                        double *density_int_cols,
                                                        int     num_layers,
                                                        int     num_columns) {
    int layer_idx = blockIdx.x * blockDim.x + threadIdx.x;
    //    int col_idx   = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;

    if (layer_idx < num_layers) {
        // get offset into start of column data
        double *delta_col_mass_upper = &(delta_col_mass_upper_cols[col_block_idx * num_layers]);
        double *delta_col_mass_lower = &(delta_col_mass_lower_cols[col_block_idx * num_layers]);

        double *density_int = &(density_int_cols[col_block_idx * (num_layers + 1)]);
        double *density_lay = &(density_lay_cols[col_block_idx * num_layers]);

        delta_col_mass_upper[layer_idx] = 0.5
                                          * (density_int[layer_idx + 1] + density_lay[layer_idx])
                                          * (altitude_int[layer_idx + 1] - altitude_lay[layer_idx]);

        delta_col_mass_lower[layer_idx] = 0.5 * (density_int[layer_idx] + density_lay[layer_idx])
                                          * (altitude_lay[layer_idx] - altitude_int[layer_idx]);
    }
}

// initialise delta_colmass arrays from pressure
// same as helios.source.host_functions.construct_grid
__global__ void initialise_delta_colmass_density_iso(double *delta_col_mass_cols,
                                                     double *altitude_int,
                                                     double *density_lay_cols,
                                                     int     num_layers,
                                                     int     num_columns) {
    int layer_idx = blockIdx.x * blockDim.x + threadIdx.x;
    // int col_idx   = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;
    int col_size      = num_layers;

    if (layer_idx < num_layers) {
        // get offset into start of column data
        double *delta_col_mass = &(delta_col_mass_cols[col_block_idx * col_size]);

        double *density_lay = &(density_lay_cols[col_block_idx * num_layers]);
        delta_col_mass[layer_idx] =
            density_lay[layer_idx] * (altitude_int[layer_idx + 1] - altitude_int[layer_idx]);
    }
}

// ************************************************************************************************

// single column pressure and temperature interpolation from layers to
// interfaces needs to loop from 0 to number of interfaces (nvi = nv+1) same as
// profX_RT
__global__ void interpolate_temperature_and_pressure(double *temperature_lay_cols,      // out
                                                     double *temperature_lay_thor_cols, // in
                                                     double *temperature_int_cols,      // out
                                                     double *pressure_lay_cols,         // in
                                                     double *pressure_int_cols,         // out
                                                     double *density_lay_cols,          // in
                                                     double *density_int_cols,          // out
                                                     double *altitude_lay,              // in
                                                     double *altitude_int,              // in
                                                     double *Tsurface_cols,
                                                     double  T_intern,
                                                     double  gravit,
                                                     int     num_layers,
                                                     int     num_columns) {
    int int_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int col_idx = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;

    if (col_idx < num_columns) {
        // Get offset arrays into columns
        double *temperature_lay      = &(temperature_lay_cols[col_block_idx * (num_layers + 1)]);
        double *temperature_lay_thor = &(temperature_lay_thor_cols[col_block_idx * num_layers]);
        double *temperature_int      = &(temperature_int_cols[col_block_idx * (num_layers + 1)]);
        double *pressure_lay         = &(pressure_lay_cols[col_block_idx * num_layers]);
        double *pressure_int         = &(pressure_int_cols[col_block_idx * (num_layers + 1)]);
        double *density_lay          = &(density_lay_cols[col_block_idx * num_layers]);
        double *density_int          = &(density_int_cols[col_block_idx * (num_layers + 1)]);
        double *Tsurface             = &(Tsurface_cols[col_block_idx]);

        // Prepare temperature array with T_intern
        if (int_idx < num_layers) {
            temperature_lay[int_idx] = temperature_lay_thor[int_idx];
        }
        else if (int_idx == num_layers) {
            temperature_lay[num_layers] = Tsurface[0];
        }

        // compute interface values
        if (int_idx == 0) {
            // extrapolate to lower boundary
            double psm = pressure_lay[1]
                         - density_lay[0] * gravit
                               * (2 * altitude_int[0] - altitude_lay[0] - altitude_lay[1]);

            double ps = 0.5 * (pressure_lay[0] + psm);

            double dbottom = density_lay[0]
                             + (density_lay[1] - density_lay[0])
                                   / (altitude_lay[1] - altitude_lay[0])
                                   * (altitude_int[0] - altitude_lay[0]);
            if (dbottom < 0.0)
                dbottom = 0.0; // prevents pressure at the top from becoming negative

            density_int[0]  = dbottom;
            pressure_int[0] = ps;
            temperature_int[0] =
                temperature_lay[0]; // T_intern; TEST: try with isothermal lower layer
                                    // // fixes weird shifted bottom point, doesn't
                                    // change upward flux being different
        }
        else if (int_idx == num_layers) {
            // extrapolate to top boundary
            double pp = pressure_lay[num_layers - 2]
                        + (pressure_lay[num_layers - 1] - pressure_lay[num_layers - 2])
                              / (altitude_lay[num_layers - 1] - altitude_lay[num_layers - 2])
                              * (2 * altitude_int[num_layers] - altitude_lay[num_layers - 1]
                                 - altitude_lay[num_layers - 2]);
            if (pp < 0.0)
                pp = 0.0; // prevents pressure at the top from becoming negative
            double ptop = 0.5 * (pressure_lay[num_layers - 1] + pp);

            pressure_int[num_layers] = ptop;

            double dtop = density_lay[num_layers - 1]
                          + (density_lay[num_layers - 1] - density_lay[num_layers - 2])
                                / (altitude_lay[num_layers - 1] - altitude_lay[num_layers - 2])
                                * (altitude_int[num_layers] - altitude_lay[num_layers - 1]);
            if (dtop < 0.0)
                dtop = 0.0; // prevents pressure at the top from becoming negative

            density_int[num_layers] = dtop;

            // extrapolate to top interface
            temperature_int[num_layers] = temperature_lay_thor[num_layers - 1]
                                          + 0.5
                                                * (temperature_lay_thor[num_layers - 1]
                                                   - temperature_lay_thor[num_layers - 2]);
        }
        else if (int_idx < num_layers) {
            // interpolation between layers
            // Helios computes gy taking the middle between the layers. We can have
            // non uniform Z levels, so linear interpolation
            double xi       = altitude_int[int_idx];
            double xi_minus = altitude_lay[int_idx - 1];
            double xi_plus  = altitude_lay[int_idx];
            double a        = (xi - xi_plus) / (xi_minus - xi_plus);
            double b        = (xi - xi_minus) / (xi_plus - xi_minus);

            pressure_int[int_idx] = pressure_lay[int_idx - 1] * a + pressure_lay[int_idx] * b;

            temperature_int[int_idx] =
                temperature_lay_thor[int_idx - 1] * a + temperature_lay_thor[int_idx] * b;

            density_int[int_idx] = density_lay[int_idx - 1] * a + density_lay[int_idx] * b;
        }
    }
}

__global__ void
increment_Qheat(double *Qheat_global, double *Qheat, double scaling, int num_sample) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_sample) {
        // delta_flux/delta_z
        Qheat_global[idx] += scaling * Qheat[idx];
    }
}

__global__ void compute_column_Qheat(double *F_net_cols, // net flux, layer
                                     double *z_int,
                                     double *Qheat_cols,
                                     double *Tsurface_cols,
                                     double  F_intern,
                                     double  Csurf,
                                     double  time_step,
                                     int     num_layers,
                                     int     tot_num_cols,
                                     bool    surface) {
    int layer_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int col_idx   = blockIdx.z * blockDim.z + threadIdx.z;
    // index of column in column batch.
    int col_block_idx = blockIdx.z;
    if (col_idx < tot_num_cols) {
        double *F_net    = &(F_net_cols[col_block_idx * (num_layers + 1)]);
        double *Qheat    = &(Qheat_cols[col_idx * num_layers]);
        double *Tsurface = &(Tsurface_cols[col_block_idx]);

        if (layer_idx == 0) {
            // delta_flux/delta_z
            // F_net positive in upward direction (F_up - F_down)
            // F_intern positive, flux out of bottom surface
            // Qheat negative when net flux differential out of layer is positive
            if (surface) {
                Tsurface[0] += (-F_net[0] + F_intern) * time_step / Csurf;
                Qheat[layer_idx] = -((F_net[1] - (F_net[0]))) / (z_int[1] - z_int[0]);
            }
            else {
                Qheat[layer_idx] = -((F_net[1] - (F_net[0] + F_intern))) / (z_int[1] - z_int[0]);
            }
        }
        else if (layer_idx < num_layers) {
            // delta_flux/delta_z
            Qheat[layer_idx] = -(F_net[layer_idx + 1] - F_net[layer_idx])
                               / (z_int[layer_idx + 1] - z_int[layer_idx]);
        }
    }
}

bool two_streams_radiative_transfer::phy_loop(ESP &                  esp,
                                              const SimulationSetup &sim,
                                              kernel_diagnostics &   diag,
                                              int                    nstep, // Step number
                                              double                 time_step)             // Time-step [s]
{
    bool run      = true;
    qheat_scaling = 1.0;

    if (spinup_start_step > -1 && spinup_stop_step > -1) {
        if (nstep < spinup_start_step) // before spinup
        {
            run           = false;
            qheat_scaling = 0.0;
        }
        else if ((nstep >= spinup_start_step) && (nstep <= spinup_stop_step)) // during spinup
        {
            double x = (double)(nstep - spinup_start_step)
                       / (double)(spinup_stop_step - spinup_start_step);
            qheat_scaling = (1 + sin(M_PI * x - M_PI / 2.0)) / 2.0;
            run           = true;
        }
    }

    if (spindown_start_step > -1 && spindown_stop_step > -1) {
        if ((nstep >= spindown_start_step) && (nstep <= spindown_stop_step)) {
            double x = (double)(nstep - spindown_start_step)
                       / (double)(spindown_stop_step - spindown_start_step);
            qheat_scaling = 1.0 - (1 + sin(M_PI * x - M_PI / 2.0)) / 2.0;
            run           = true;
        }
        else if (nstep >= spindown_stop_step) {
            run           = false;
            qheat_scaling = 0.0;
        }
    }

    if (run) {

        alf.debug_nstep = nstep;

        const int num_blocks = 256;
        const int num_cols   = num_parallel_columns;
        if (nstep % compute_every_n_iteration == 0 || start_up) {
            std::shared_ptr<double[]> col_cos_zenith_angle_h =
                esp.insolation.get_host_cos_zenith_angles();
            double *col_cos_zenith_angle_d = esp.insolation.get_device_cos_zenith_angles();
            Qheat.zero();
            F_down_tot.zero();
            F_up_tot.zero();
            F_dir_tot.zero();
            F_up_band.zero();
            F_dir_band.zero();
            F_net.zero();

            g0_tot.zero();
            w0_tot.zero();

            printf("\r\n");
            printf("\r\n");
            printf("\r\n");
            cudaDeviceSynchronize();
            cuda_check_status_or_exit(__FILE__, __LINE__);
            int nbin = alf.opacities.nbin;
            // loop on columns
            for (int column_idx = 0; column_idx < esp.point_num;
                 column_idx += num_parallel_columns) {
                // printf("column_idx_tsrt: %d\n", column_idx);

                int current_num_cols = min(num_cols, esp.point_num - column_idx);
                alf.debug_col_idx    = column_idx;

#ifdef COLUMN_LOOP_PROGRESS_BAR
                print_progress((column_idx + 1.0) / double(esp.point_num));
#endif // COLUMN_LOOP_PROGRESS_BAR

                F_up_wg.zero();
                F_down_wg.zero();
                F_dir_wg.zero();
                if (iso) {
                }
                else {
                    if (!thomas) {
                        Fc_down_wg.zero();
                        Fc_up_wg.zero();
                    }
                    Fc_dir_wg.zero();
                }

                alf.reset();

                pressure_int.zero();
                temperature_int.zero();
                temperature_lay.zero();
                density_int.zero();

                int num_layers = esp.nv;

                int column_offset = column_idx * num_layers;

                double gravit = sim.Gravit;
                // fetch column values

                double *column_layer_temperature_thor = &(esp.temperature_d[column_offset]);
                double *column_layer_pressure         = &(esp.pressure_d[column_offset]);
                double *column_density                = &(esp.Rho_d[column_offset]);
                // initialise interpolated T and P
                double *Tsurface_col =
                    &(esp.Tsurface_d[column_idx]); //when !surface, this array is empty

                // use mu_star per column

#ifdef DUMP_HELIOS_TP
                cudaDeviceSynchronize();
                cuda_check_status_or_exit(__FILE__, __LINE__);

                {
                    // get data for current column batch TP profile
                    std::shared_ptr<double[]> pressure_h =
                        get_cuda_data(column_layer_pressure, esp.nv * current_num_cols);
                    std::shared_ptr<double[]> temperature_h =
                        get_cuda_data(column_layer_temperature_thor, esp.nv * current_num_cols);
                    std::shared_ptr<double[]> density_h =
                        get_cuda_data(column_density, esp.nv * current_num_cols);

                    for (int c = 0; c < current_num_cols; c++) {
                        // dump a TP profile for HELIOS input
                        if ((column_idx + c) % HELIOS_TP_STRIDE == 0) {
                            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                                         + "/alfprof/"
                                                           "step_"
                                                         + std::to_string(nstep) + "/column_"
                                                         + std::to_string(column_idx + c) + "/";
                            create_output_dir(DBG_OUTPUT_DIR);
                            double lon = esp.lonlat_h[(column_idx + c) * 2 + 0] * 180 / M_PI;
                            double lat = esp.lonlat_h[(column_idx + c) * 2 + 1] * 180 / M_PI;

                            double p_toa   = pressure_h[esp.nv * c + esp.nv - 1];
                            double p_boa   = pressure_h[esp.nv * c];
                            double mu_star = -col_cos_zenith_angle_h[column_idx + c];

                            // Print out initial TP profile
                            string output_file_name = DBG_OUTPUT_DIR + "tpprofile_init.dat";

                            FILE * tp_output_file = fopen(output_file_name.c_str(), "w");
                            string comment =
                                "# Helios TP profile table at lat: [" + std::to_string(lon)
                                + "] lon: [" + std::to_string(lat) + "] mustar: ["
                                + std::to_string(mu_star) + "] P_BOA: [" + std::to_string(p_boa)
                                + "] P_TOA: [" + std::to_string(p_toa) + "]\n";

                            fprintf(tp_output_file, comment.c_str());
                            fprintf(tp_output_file, "#\tT[K]\tP[bar]\trho\n");

                            for (int i = 0; i < esp.nv; i++) {
                                fprintf(tp_output_file,
                                        "%#.6g\t%#.6g\t%#.6g\n",
                                        temperature_h[i + esp.nv * c],
                                        pressure_h[i + esp.nv * c] / 1e5,
                                        density_h[i + esp.nv * c]);
                            }

                            fclose(tp_output_file);
                        }
                    }
                }
#endif // DUMP_HELIOS_TP
                cudaDeviceSynchronize();
                cuda_check_status_or_exit(__FILE__, __LINE__);
                {
                    dim3 grid(int((num_layers + 1) / num_blocks) + 1, 1, num_cols);
                    dim3 block(num_blocks, 1, 1);
                    interpolate_temperature_and_pressure<<<grid, block>>>(
                        *temperature_lay,
                        column_layer_temperature_thor,
                        *temperature_int,
                        column_layer_pressure,
                        *pressure_int,
                        column_density,
                        *density_int,
                        esp.Altitude_d,
                        esp.Altitudeh_d,
                        Tsurface_col,
                        T_internal,
                        gravit,
                        num_layers,
                        current_num_cols);
                    cudaDeviceSynchronize();
                    cuda_check_status_or_exit(__FILE__, __LINE__);
                }
                BENCH_POINT_I_S(
                    nstep, column_idx, "Alf_interpTnP", (), ("T_lay", "T_int", "P_int"));

#ifdef DUMP_HELIOS_TP
                {
                    // dump a TP profile for HELIOS input
                    std::shared_ptr<double[]> pressure_int_h    = pressure_int.get_host_data();
                    std::shared_ptr<double[]> temperature_int_h = temperature_int.get_host_data();
                    std::shared_ptr<double[]> density_int_h     = density_int.get_host_data();

                    for (int c = 0; c < current_num_cols; c++) {

                        if (column_idx % HELIOS_TP_STRIDE == 0) {
                            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                                         + "/alfprof/"
                                                           "step_"
                                                         + std::to_string(nstep) + "/column_"
                                                         + std::to_string(column_idx + c) + "/";
                            create_output_dir(DBG_OUTPUT_DIR);

                            double lon = esp.lonlat_h[(column_idx + c) * 2 + 0] * 180 / M_PI;
                            double lat = esp.lonlat_h[(column_idx + c) * 2 + 1] * 180 / M_PI;

                            // get col mu star from zenith angle

                            double p_toa = pressure_int_h[esp.nvi * c + esp.nvi - 1];
                            double p_boa = pressure_int_h[esp.nvi * c];

                            double mu_star = -col_cos_zenith_angle_h[column_idx + c];
                            // Print out initial TP profile
                            string output_file_name = DBG_OUTPUT_DIR + "tpprofile_interface.dat";

                            FILE * tp_output_file = fopen(output_file_name.c_str(), "w");
                            string comment        = "# Helios TP interface profile table at lat: ["
                                             + std::to_string(lon) + "] lon: ["
                                             + std::to_string(lat) + "] mustar: ["
                                             + std::to_string(mu_star) + "] P_BOA: ["
                                             + std::to_string(p_boa) + "] P_TOA: ["
                                             + std::to_string(p_toa) + "]\n";

                            fprintf(tp_output_file, comment.c_str());
                            fprintf(tp_output_file, "#\tT[K]\tP[bar]\trho\n");

                            for (int i = 0; i < esp.nvi; i++) {
                                fprintf(tp_output_file,
                                        "%#.6g\t%#.6g\t%#.6g\n",
                                        temperature_int_h[i + esp.nvi * c],
                                        pressure_int_h[i + esp.nvi * c] / 1e5,
                                        density_int_h[i + esp.nvi * c]);
                            }

                            fclose(tp_output_file);
                        }
                    }
                }
#endif // DUMP_HELIOS_TP

                // initialise delta_col_mass
                if (false) { //this is not exectuted (old hydrostatic def'n of column mass)
                    // initialise delta_col_mass
                    if (iso) {
                        dim3 grid(int((num_layers + 1) / num_blocks) + 1, 1, current_num_cols);
                        dim3 block(num_blocks, 1, 1);
                        initialise_delta_colmass_pressure_iso<<<grid, block>>>(*alf.delta_col_mass,
                                                                               *pressure_int,
                                                                               gravit,
                                                                               num_layers,
                                                                               current_num_cols);
                    }
                    else {
                        dim3 grid(int((num_layers + 1) / num_blocks) + 1, 1, current_num_cols);
                        dim3 block(num_blocks, 1, 1);
                        initialise_delta_colmass_pressure_noniso<<<grid, block>>>(
                            *alf.delta_col_upper,
                            *alf.delta_col_lower,
                            column_layer_pressure,
                            *pressure_int,
                            gravit,
                            num_layers,
                            current_num_cols);
                    }
                }
                else { //non-hydrostatic definition of column mass
                    if (iso) {
                        dim3 grid(int((num_layers + 1) / num_blocks) + 1, 1, current_num_cols);
                        dim3 block(num_blocks, 1, 1);
                        initialise_delta_colmass_density_iso<<<grid, block>>>(*alf.delta_col_mass,
                                                                              esp.Altitudeh_d,
                                                                              column_density,
                                                                              num_layers,
                                                                              current_num_cols);
                    }
                    else {
                        dim3 grid(int((num_layers + 1) / num_blocks) + 1, 1, current_num_cols);
                        dim3 block(num_blocks, 1, 1);
                        initialise_delta_colmass_density_noniso<<<grid, block>>>(
                            *alf.delta_col_upper,
                            *alf.delta_col_lower,
                            esp.Altitude_d,
                            esp.Altitudeh_d,
                            column_density,
                            *density_int,
                            num_layers,
                            current_num_cols);
                    }
                }
                cudaDeviceSynchronize();
                cuda_check_status_or_exit(__FILE__, __LINE__);

                double *z_lay = esp.Altitude_d;
                double *z_int = esp.Altitudeh_d;
                // internal to alfrodull_engine

                double *dev_starflux = *star_flux;
                // limit where to switch from noniso to iso equations to keep model
                // stable as defined in host_functions.set_up_numerical_parameters
                double delta_tau_limit = 1e-4;

                // compute fluxes

                // singlewalk (for iterative solver)
                //  true -> 201 iterations,
                //  false -> 4 iterations,

                bool singlewalk_loc = scat_single_walk;
                int  ninterface     = nlayer + 1;

                {
                    double *cos_zenith_angle_cols = &(col_cos_zenith_angle_d[column_idx]);
                    int     column_offset_int     = column_idx * ninterface;

                    double *F_col_down_tot = &((*F_down_tot)[column_offset_int]);
                    double *F_col_up_tot   = &((*F_up_tot)[column_offset_int]);
                    double *F_col_dir_tot  = &((*F_dir_tot)[column_offset_int]);
                    double *F_col_net      = &((*F_net)[column_offset_int]);

                    double *F_col_dir_band = nullptr;
                    if (store_dir_spectrum)
                        F_col_dir_band = &((*F_dir_band)[column_idx * ninterface * nbin]);
                    else
                        F_col_dir_band = *F_dir_band;
                    double *F_up_TOA_spectrum_col = &((*F_up_TOA_spectrum)[column_idx * nbin]);

                    alf.compute_radiative_transfer(dev_starflux,          // dev_starflux
                                                   *temperature_lay,      // dev_T_lay
                                                   *temperature_int,      // dev_T_int
                                                   column_layer_pressure, // dev_p_lay
                                                   *pressure_int,         // dev_p_int
                                                   false,                 // interp_press_and_temp
                                                   true,           // interp_and_calc_flux_step
                                                   z_lay,          // z_lay
                                                   singlewalk_loc, // singlewalk
                                                   *F_down_wg,
                                                   *F_up_wg,
                                                   *Fc_down_wg,
                                                   *Fc_up_wg,
                                                   *F_dir_wg,
                                                   *Fc_dir_wg,
                                                   delta_tau_limit,
                                                   F_col_down_tot,
                                                   F_col_up_tot,
                                                   F_col_dir_tot,
                                                   F_col_net,
                                                   *F_down_band,
                                                   *F_up_band,
                                                   F_col_dir_band,
                                                   F_up_TOA_spectrum_col,
                                                   cos_zenith_angle_cols,
                                                   *surface_albedo,
                                                   current_num_cols,
                                                   column_idx,
                                                   esp.surface);
                    cudaDeviceSynchronize();
                    cuda_check_status_or_exit(__FILE__, __LINE__);
                }

                // get the g0 and w0 integrated
                if (store_w0_g0) {
                    double *g0_tot_col = &((*g0_tot)[column_idx * nlayer * nbin]);
                    double *w0_tot_col = &((*w0_tot)[column_idx * nlayer * nbin]);
                    alf.get_column_integrated_g0_w0(g0_tot_col, w0_tot_col, current_num_cols);
                }

                // compute Delta flux

                {
                    dim3    grid(int((esp.nv / num_blocks) + 1), 1, current_num_cols);
                    dim3    block(num_blocks, 1, 1);
                    double *qheat     = &((*Qheat)[column_idx * nlayer]);
                    double *F_col_net = &((*F_net)[column_idx * ninterface]);
                    compute_column_Qheat<<<grid, block>>>(F_col_net, // net flux, layer
                                                          z_int,
                                                          qheat,
                                                          Tsurface_col,
                                                          F_intern,
                                                          esp.Csurf,
                                                          time_step,
                                                          num_layers,
                                                          esp.point_num,
                                                          esp.surface);
                    cudaDeviceSynchronize();
                    cuda_check_status_or_exit(__FILE__, __LINE__);
                }
#ifdef DEBUG_PRINTOUT_ARRAYS
                debug_print_columns(
                    esp, col_cos_zenith_angle_h, nstep, column_idx, current_num_cols);
#endif // DEBUG_PRINTOUT_ARRAYS
            }
            start_up = false;
        }

        printf("\r\n");

        cudaDeviceSynchronize();
        cuda_check_status_or_exit(__FILE__, __LINE__);
        // set Qheat
        int num_samples = (esp.point_num * nlayer);
        increment_Qheat<<<(num_samples / num_blocks) + 1, num_blocks>>>(
            esp.profx_Qheat_d, *Qheat, qheat_scaling, num_samples);
        cudaDeviceSynchronize();
        cuda_check_status_or_exit(__FILE__, __LINE__);
    }
    last_step = nstep;

    BENCH_POINT_I(nstep, "Alf_phy_loop_E", (), ("F_up_tot", "F_down_tot", "AlfQheat"));

    return true;
}

bool two_streams_radiative_transfer::store_init(storage &s) {
    if (!s.has_table("/Tstar"))
        s.append_value(T_star, "/Tstar", "K", "Temperature of host star");
    // s.append_value(Tint, "/alf_Tint", "K", "Temperature of interior heat
    // flux");
    if (!s.has_table("/planet_star_dist"))
        s.append_value(planet_star_dist_config,
                       "/planet_star_dist",
                       "au",
                       "distance b/w host star and planet");

    if (!s.has_table("/radius_star"))
        s.append_value(R_star_config, "/radius_star", "R_sun", "radius of host star");

    s.append_value(thomas ? 1 : 0, "/alf_thomas", "-", "Alfrodull use Thomas Solver");
    s.append_value(iso ? 1.0 : 0.0, "/alf_isothermal", "-", "Isothermal layers");
    s.append_value(
        real_star ? 1.0 : 0.0, "/alf_real_star", "-", "Alfrodull use real star spectrum or Planck");

    s.append_value(
        fake_opac ? 1.0 : 0.0, "/alf_fake_opac", "-", "Alfrodull use artificial opacity");
    s.append_value(scat ? 1.0 : 0.0, "/alf_scat", "-", "Alfrodull scattering");

    s.append_value(scat_corr ? 1.0 : 0.0,
                   "/alf_scat_corr",
                   "-",
                   "Alfrodull Improved two-stream scattering correction");

    s.append_value(scat_single_walk ? 1 : 0,
                   "/alf_scat_single_walk",
                   "-",
                   "Iterative solver single walk mode");
    s.append_value(g_0, "/alf_g_0", "-", "asymmetry factor");
    s.append_value(diffusivity, "/alf_diffusivity", "-", "Diffusivity factor");
    s.append_value(epsi, "/alf_epsilon", "-", "One over Diffusivity factor");
    s.append_value(epsilon_2, "/alf_epsilon_2", "-", "Epsilon 2 factor");

    s.append_value(experimental_opacities_offset,
                   "/Alf_experimental_opacities_offset",
                   "-",
                   "Alfrodull Experimental opacity offset for debugging");

    s.append_value(alf.opacities.ny, "/alf_ny", "-", "Alfrodull number of weights in bins");

    s.append_value(i2s_transition, "/alf_i2s_transition", "-", "Alfrodull i2s transition");

    // s.append_value(G_pm_limiter ? 1 : 0,
    //                "/G_pm_max_limiter",
    //                "-",
    //                "Alfrodull limiter on abs(G_pm) at 1e8 (HELIOS compatible)");
    // s.append_value(G_pm_denom_limit,
    //                "/G_pm_factor_limit",
    //                "-",
    //                "Alfrodull limiter on abs(G_pm) or G_pm_denominator");
    //
    // s.append_value(G_pm_limit_on_full_G_pm ? 1 : 0,
    //                "/G_pm_limit_on_full_G_pm",
    //                "-",
    //                "Alf G_pm limit applies to full G_pm term or only to G_pm denominator");
    //
    // s.append_value(mu_star_wiggle_increment,
    //                "/Alf_G_pm_mu_star_increment",
    //                "-",
    //                "Alf increment in degrees applied to incomming angle when hit "
    //                "G_pm limit");
    // s.append_value(wiggle_iteration_max,
    //                "/Alf_mu_star_iteration_max",
    //                "-",
    //                "Alf max number of incoming angle increases on G_pm limit");
    s.append_value(mu_star_limit_degrees,
                   "/Alf_direct_beam_angle_limit",
                   "-",
                   "Alf limit on incoming angle to tangent to set direct beam to 0.0");
    // commented out, no easy way to store strings in the table in the same way.
    // if (real_star) {
    //     std::string str[] = {stellar_spectrum_file};
    //     s.append_value(str, "/Alf_stellar_spectrum", "-", "spectrum file
    //     name");
    // }
    // {
    //     std::string str[] = {opacities_file};
    //     s.append_value(str, "/Alf_opacities_file", "-", "Opacities file name");
    // }
    s.append_value(
        opacity_file_is_CGS ? 1 : 0, "/Alf_opacities_file_in_CGS", "-", "Opacities input in CGS");
    s.append_value(spinup_start_step, "/Alf_spinup_start", "-", "Alf spinup start step");
    s.append_value(spinup_stop_step, "/Alf_spinup_stop", "-", "Alf spinup stop step");
    s.append_value(spindown_start_step, "/Alf_spindown_start", "-", "Alf spindown start step");
    s.append_value(spindown_stop_step, "/Alf_spindown_stop", "-", "Alf spindown stop step");

    s.append_value(num_parallel_columns,
                   "/Alf_num_parallel_columns",
                   "-",
                   " Alf number of columns to compute in parallel");
    s.append_value(debug_output ? 1 : 0, "/Alf_debug", "-", "Alf Debug output");
    s.append_value(store_w0_g0 ? 1 : 0, "/Alf_store_w0_g0", "-", "Alf store w0 g0 per band");
    s.append_value(store_dir_spectrum ? 1 : 0,
                   "/Alf_store_dir_spectrum",
                   "-",
                   "Alf store directional beam per band");
    s.append_value(null_planck_function ? 1 : 0,
                   "/Alf_null_planck_function",
                   "-",
                   "Alf set Planck function to 0");
    // {
    //     std::string str[] = {cloud_filename};
    //     s.append_value(str, "/Alf_cloudfile", "-", "Alf cloud opacity file");
    // }
    s.append_value(compute_every_n_iteration,
                   "/alf_compute_periodicity",
                   "n",
                   "Alfrodull compute periodicity");
    // s.append_value(opacities_file, "/alf_opacity_file", "path", "Alfrodull
    // opacitiy file used");

    s.append_value(dir_beam ? 1.0 : 0.0, "/alf_dir_beam", "-", "Direct irradiation beam");
    s.append_value(geom_zenith_corr ? 1.0 : 0.0,
                   "/alf_geom_zenith_corr",
                   "-",
                   "Geometric zenith angle correction");

    s.append_value(alf.opacities.nbin, "/alf_num_bands", "-", "Number of wavelength_bands for Alf");
    s.append_value(
        store_w0_g0 ? 1 : 0, "/alf_w0_g0_per_band", "-", "Stored w0 and g0 per band for Alf");

    s.append_value(clouds ? 1 : 0, "/alf_cloud", "-", "Simulate clouds");
    s.append_value(fcloud, "/alf_fcloud", "-", "f_cloud");

    return true;
}
//***************************************************************************************************

bool two_streams_radiative_transfer::store(const ESP &esp, storage &s) {
    std::shared_ptr<double[]> F_net_h = F_net.get_host_data();
    s.append_table(F_net_h.get(), F_net.get_size(), "/F_net", "W m^-2", "Net Flux");

    std::shared_ptr<double[]> Qheat_h = Qheat.get_host_data();
    s.append_table(Qheat_h.get(), Qheat.get_size(), "/Alf_Qheat", "W m^-3", "Alfrodull Qheat");

    std::shared_ptr<double[]> F_up_tot_h = F_up_tot.get_host_data();
    s.append_table(
        F_up_tot_h.get(), F_up_tot.get_size(), "/F_up_tot", "W m^-2", "Total upward flux");

    std::shared_ptr<double[]> F_down_tot_h = F_down_tot.get_host_data();
    s.append_table(
        F_down_tot_h.get(), F_down_tot.get_size(), "/F_down_tot", "W m^-2", "Total downward flux");

    std::shared_ptr<double[]> F_dir_tot_h = F_dir_tot.get_host_data();
    s.append_table(
        F_dir_tot_h.get(), F_dir_tot.get_size(), "/F_dir_tot", "W m^-2", "Total beam flux");

    if (store_dir_spectrum) {
        std::shared_ptr<double[]> F_dir_band_h = F_dir_band.get_host_data();
        s.append_table(F_dir_band_h.get(),
                       F_dir_band.get_size(),
                       "/F_dir_band",
                       "W m^-2 m^-1",
                       "Directional beam spectrum");
    }

    if (store_w0_g0) {
        std::shared_ptr<double[]> w0_tot_h = w0_tot.get_host_data();
        s.append_table(w0_tot_h.get(),
                       w0_tot.get_size(),
                       "/w0_band",
                       " ",
                       "Single scattering albedo per band");

        std::shared_ptr<double[]> g0_tot_h = g0_tot.get_host_data();
        s.append_table(g0_tot_h.get(), g0_tot.get_size(), "/g0_band", " ", "asymmetry per band");
    }

    {
        int                       nbin             = alf.opacities.nbin;
        int                       numinterfaces    = esp.nvi;
        std::shared_ptr<double[]> planckband_lay_h = alf.planckband_lay.get_host_data();
        std::shared_ptr<double[]> spectrum         = std::shared_ptr<double[]>(new double[nbin]);
        for (int i = 0; i < nbin; i++)
            spectrum[i] = planckband_lay_h[(numinterfaces - 1) + i * (numinterfaces - 1 + 2)];

        s.append_table(spectrum.get(),
                       nbin,
                       "/alf_stellar_spectrum",
                       "W m^-2 m^-1",
                       "Alfrodull stellar spectrum");
    }

    std::shared_ptr<double[]> F_up_TOA_spectrum_h = F_up_TOA_spectrum.get_host_data();
    s.append_table(F_up_TOA_spectrum_h.get(),
                   F_up_TOA_spectrum.get_size(),
                   "/F_up_TOA_spectrum",
                   "W m^-2 m^-1",
                   "Upward Flux per bin at TOA");

    std::shared_ptr<double[]> lambda_wave_h = alf.opacities.dev_opac_wave.get_host_data();
    s.append_table(lambda_wave_h.get(),
                   alf.opacities.dev_opac_wave.get_size(),
                   "/lambda_wave",
                   "m",
                   "Center wavelength");

    std::shared_ptr<double[]> lambda_interwave_h = alf.opacities.dev_opac_interwave.get_host_data();
    s.append_table(lambda_interwave_h.get(),
                   alf.opacities.dev_opac_interwave.get_size(),
                   "/lambda_interwave",
                   "m",
                   "Interface wavelength");

    std::shared_ptr<double[]> lambda_deltawave_h = alf.opacities.dev_opac_deltawave.get_host_data();
    s.append_table(lambda_deltawave_h.get(),
                   alf.opacities.dev_opac_deltawave.get_size(),
                   "/lambda_deltawave",
                   "m",
                   "Wavelength width of bins");

    s.append_value(qheat_scaling, "/qheat_scaling", "-", "QHeat scaling");

    return true;
}

bool two_streams_radiative_transfer::free_memory() {
    return true;
}

// ***************************************************************************************************************
// ***************************************************************************************************************
// special version for Thomas algorithm X_buff that has special indexing
// works only for ny = 1
void two_streams_radiative_transfer::print_X_buff_thomas_data_to_file(
    ESP &                       esp,
    int                         nstep,
    int                         column_idx,
    int                         num_stack,
    int                         num_cols,
    string                      stackname,
    cuda_device_memory<double> &array_,
    string                      output_file_base) {

    double *array = *array_;

    int array_size = array_.get_size();
    int nbin       = alf.opacities.nbin;
    int ny         = alf.opacities.ny;

    // Print out single scattering albedo data

    if (ny != 1)
        printf("Debug print of X_thomas error, ny != 1\n");

    // std::shared_ptr<double[]> array_h =
    //     integrate_band(array, *alf.gauss_weights, num_val, nbin, ny);
    std::shared_ptr<double[]> array_h = get_cuda_data(array, array_size);

    cuda_check_status_or_exit((string(__FILE__ ":") + string(output_file_base)).c_str(), __LINE__);

    std::shared_ptr<double[]> delta_lambda_h =
        get_cuda_data(*alf.opacities.dev_opac_deltawave, nbin);

    cuda_check_status_or_exit((string(__FILE__ ":") + string(output_file_base)).c_str(), __LINE__);

    for (int c = 0; c < num_cols; c++) {
        // offset of column
        int col_offset = c * (2 * nlayer + 1) * nbin * 2;

        // loop on columns
        // printf("/alfprof/step_%d/column_%d/\n", nstep, column_idx + c);
        std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                     + "/alfprof/"
                                       "step_"
                                     + std::to_string(nstep) + "/column_"
                                     + std::to_string(column_idx + c) + "/";
        create_output_dir(DBG_OUTPUT_DIR);

        string output_file_name = DBG_OUTPUT_DIR + output_file_base + ".dat";

        FILE *output_file = fopen(output_file_name.c_str(), "w");
        fprintf(output_file, "bin\t");
        fprintf(output_file, "deltalambda\t");
        for (int i = 0; i < num_stack; i++)
            fprintf(output_file, "%s[%d]\t", stackname.c_str(), i);
        fprintf(output_file, "\n");

        for (int b = 0; b < nbin; b++) {
            fprintf(output_file, "%d\t", b);
            fprintf(output_file, "%#.6g\t", delta_lambda_h[b]);
            for (int i = 0; i < num_stack; i++) {
                fprintf(output_file, "%#.6g\t", array_h[i + col_offset + (2 * nlayer + 1) * b * 2]);
            }
            fprintf(output_file, "\n");
        }
        fclose(output_file);
    }
}
// ***************************************************************************************************************
void two_streams_radiative_transfer::print_weighted_band_data_to_file(
    ESP &                       esp,
    int                         nstep,
    int                         column_idx,
    int                         num_stack,
    int                         num_cols,
    string                      stackname,
    cuda_device_memory<double> &array,
    string                      output_file_base,
    bool                        global) {
    print_weighted_band_data_to_file(esp,
                                     nstep,
                                     column_idx,
                                     num_stack,
                                     num_cols,
                                     stackname,
                                     *array,
                                     array.get_size(),
                                     output_file_base,
                                     global);
}

void two_streams_radiative_transfer::print_weighted_band_data_to_file(ESP &   esp,
                                                                      int     nstep,
                                                                      int     column_idx,
                                                                      int     num_stack,
                                                                      int     num_cols,
                                                                      string  stackname,
                                                                      double *array,
                                                                      int     array_size,
                                                                      string  output_file_base,
                                                                      bool    global) {
    int nbin = alf.opacities.nbin;
    int ny   = alf.opacities.ny;

    // Print out single scattering albedo data

    int num_val = array_size / (nbin * ny);

    std::shared_ptr<double[]> array_h =
        integrate_band(array, *alf.gauss_weights, num_val, nbin, ny);
    cuda_check_status_or_exit((string(__FILE__ ":") + string(output_file_base)).c_str(), __LINE__);

    for (int c = 0; c < num_cols; c++) {
        int access_idx = global ? c + column_idx : c;

        // loop on columns
        std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                     + "/alfprof/"
                                       "step_"
                                     + std::to_string(nstep) + "/column_"
                                     + std::to_string(column_idx + c) + "/";
        create_output_dir(DBG_OUTPUT_DIR);

        string output_file_name = DBG_OUTPUT_DIR + output_file_base + ".dat";

        FILE *output_file = fopen(output_file_name.c_str(), "w");

        std::shared_ptr<double[]> delta_lambda_h =
            get_cuda_data(*alf.opacities.dev_opac_deltawave, nbin);
        cuda_check_status_or_exit((string(__FILE__ ":") + string(output_file_base)).c_str(),
                                  __LINE__);
        fprintf(output_file, "bin\t");
        fprintf(output_file, "deltalambda\t");
        for (int i = 0; i < num_stack; i++)
            fprintf(output_file, "%s[%d]\t", stackname.c_str(), i);
        fprintf(output_file, "\n");

        for (int b = 0; b < nbin; b++) {
            fprintf(output_file, "%d\t", b);
            fprintf(output_file, "%#.6g\t", delta_lambda_h[b]);
            for (int i = 0; i < num_stack; i++) {
                fprintf(
                    output_file, "%#.6g\t", array_h[b + i * nbin + access_idx * nbin * num_stack]);
            }
            fprintf(output_file, "\n");
        }
        fclose(output_file);
    }
}

void two_streams_radiative_transfer::print_data_to_file(ESP &                       esp,
                                                        int                         nstep,
                                                        int                         column_idx,
                                                        int                         num_stack,
                                                        int                         num_cols,
                                                        string                      stackname,
                                                        string                      column_name,
                                                        cuda_device_memory<double> &array,
                                                        string output_file_base,
                                                        bool   global,
                                                        double scaling) {
    // Print out single scattering albedo data
    int                       num_val = array.get_size();
    std::shared_ptr<double[]> array_h = array.get_host_data();

    for (int c = 0; c < num_cols; c++) {
        int access_idx = global ? c + column_idx : c;
        // loop on columns
        std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                     + "/alfprof/"
                                       "step_"
                                     + std::to_string(nstep) + "/column_"
                                     + std::to_string(column_idx + c) + "/";
        create_output_dir(DBG_OUTPUT_DIR);

        string output_file_name = DBG_OUTPUT_DIR + output_file_base + ".dat";

        FILE *output_file = fopen(output_file_name.c_str(), "w");

        fprintf(output_file, "%s\t", stackname.c_str());
        fprintf(output_file, "%s\n", column_name.c_str());

        for (int i = 0; i < num_stack; i++) {
            fprintf(output_file, "%d\t%#.6g\n", i, array_h[i + access_idx * num_stack] * scaling);
        }

        fclose(output_file);
    }
    cuda_check_status_or_exit((string(__FILE__ ":") + string(output_file_base)).c_str(), __LINE__);
}

// ***************************************************************************************************************
// Helper function to print out all datasets for debugging and comparisong to
// HELIOS
void two_streams_radiative_transfer::debug_print_columns(ESP &                      esp,
                                                         std::shared_ptr<double[]> &cmustar,
                                                         int                        nstep,
                                                         int                        column_base_idx,
                                                         int                        num_cols) {
    int nbin = alf.opacities.nbin;

    {
        std::shared_ptr<double[]> temperature_h = temperature_int.get_host_data();
        std::shared_ptr<double[]> pressure_h    = pressure_int.get_host_data();

        cuda_check_status_or_exit(string(__FILE__ ":"
                                                  "tpprofile")
                                      .c_str(),
                                  __LINE__);
        for (int c = 0; c < num_cols; c++) {
            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                         + "/alfprof/"
                                           "step_"
                                         + std::to_string(nstep) + "/column_"
                                         + std::to_string(column_base_idx + c) + "/";
            create_output_dir(DBG_OUTPUT_DIR);

            double lon = esp.lonlat_h[(column_base_idx + c) * 2 + 0] * 180 / M_PI;
            double lat = esp.lonlat_h[(column_base_idx + c) * 2 + 1] * 180 / M_PI;

            // Print out initial TP profile
            string output_file_name = DBG_OUTPUT_DIR + "tprofile_interp.dat";

            FILE * tp_output_file = fopen(output_file_name.c_str(), "w");
            string comment        = "# Helios TP profile table at lat: [" + std::to_string(lon)
                             + "] lon: [" + std::to_string(lat) + "] mustar: ["
                             + std::to_string(cmustar[c]) + "]\n";

            fprintf(tp_output_file, comment.c_str());
            fprintf(tp_output_file, "#\tP(bar)\tT[K]\n");

            // fprintf(tp_output_file,
            //         "BOA\t%#.6g\n",
            //         temperature_h[esp.nv + 1 + (column_base_idx + c) * (esp.nv + 2)]);
            for (int i = 0; i <= esp.nv; i++) {
                fprintf(tp_output_file,
                        "%d\t%#.6g\t%#.6g\n",
                        i,
                        pressure_h[i + (column_base_idx + c) * (esp.nv + 2)] / 1e5,
                        temperature_h[i + (column_base_idx + c) * (esp.nv + 2)]);
            }

            fclose(tp_output_file);
        }
    }

    {
        std::shared_ptr<double[]> planck_h = alf.planckband_lay.get_host_data();

        cuda_check_status_or_exit(string(__FILE__ ":"
                                                  "plkprofile")
                                      .c_str(),
                                  __LINE__);

        for (int c = 0; c < num_cols; c++) {
            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                         + "/alfprof/"
                                           "step_"
                                         + std::to_string(nstep) + "/column_"
                                         + std::to_string(column_base_idx + c) + "/";
            create_output_dir(DBG_OUTPUT_DIR);

            // Print out planck data

            string output_file_name = DBG_OUTPUT_DIR + "plkprofile.dat";

            FILE *planck_output_file = fopen(output_file_name.c_str(), "w");

            std::shared_ptr<double[]> delta_lambda_h =
                get_cuda_data(*alf.opacities.dev_opac_deltawave, nbin);

            fprintf(planck_output_file, "bin\t");
            fprintf(planck_output_file, "deltalambda\t");
            for (int i = 0; i < esp.nv; i++)
                fprintf(planck_output_file, "layer[%d]\t", i);
            fprintf(planck_output_file, "layer[TOA]\t");
            fprintf(planck_output_file, "layer[BOA]\t");
            fprintf(planck_output_file, "\n");

            for (int b = 0; b < nbin; b++) {
                fprintf(planck_output_file, "%d\t", b);
                fprintf(planck_output_file, "%#.6g\t", delta_lambda_h[b]);
                for (int i = 0; i < esp.nv + 2; i++) {
                    fprintf(planck_output_file,
                            "%#.6g\t",
                            planck_h[b * (esp.nv + 2) + i + c * ((esp.nv + 2) * nbin)]);
                }
                fprintf(planck_output_file, "\n");
            }
            fclose(planck_output_file);
        }
    }
    if (!iso) {
        std::shared_ptr<double[]> planck_h = alf.planckband_int.get_host_data();

        cuda_check_status_or_exit(string(__FILE__ ":"
                                                  "plkintprofile")
                                      .c_str(),
                                  __LINE__);

        for (int c = 0; c < num_cols; c++) {
            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                         + "/alfprof/"
                                           "step_"
                                         + std::to_string(nstep) + "/column_"
                                         + std::to_string(column_base_idx + c) + "/";
            create_output_dir(DBG_OUTPUT_DIR);

            // Print out planck data

            string output_file_name = DBG_OUTPUT_DIR + "plkintprofile.dat";

            FILE *planck_output_file = fopen(output_file_name.c_str(), "w");

            std::shared_ptr<double[]> delta_lambda_h =
                get_cuda_data(*alf.opacities.dev_opac_deltawave, nbin);

            fprintf(planck_output_file, "bin\t");
            fprintf(planck_output_file, "deltalambda\t");
            for (int i = 0; i < esp.nvi; i++)
                fprintf(planck_output_file, "int[%d]\t", i);
            fprintf(planck_output_file, "\n");

            for (int b = 0; b < nbin; b++) {
                fprintf(planck_output_file, "%d\t", b);
                fprintf(planck_output_file, "%#.6g\t", delta_lambda_h[b]);
                for (int i = 0; i < esp.nvi; i++) {
                    fprintf(planck_output_file,
                            "%#.6g\t",
                            planck_h[b * (esp.nvi) + i + c * ((esp.nvi) * nbin)]);
                }
                fprintf(planck_output_file, "\n");
            }
            fclose(planck_output_file);
        }
    }

    {
        // Print out mean molecular weight data
        print_data_to_file(esp,
                           nstep,
                           column_base_idx,
                           esp.nv,
                           num_cols,
                           "layer",
                           "meanmolmass",
                           alf.meanmolmass_lay,
                           "meanmolmassprofile",
                           false,
                           1.0 / AMU);
    }
    if (!iso) {
        // Print out mean molecular weight data
        print_data_to_file(esp,
                           nstep,
                           column_base_idx,
                           esp.nvi,
                           num_cols,
                           "interface",
                           "meanmolmassint",
                           alf.meanmolmass_int,
                           "meanmolmassintprofile",
                           false,
                           1.0 / AMU);
    }

    if (iso) {
        print_data_to_file(esp,
                           nstep,
                           column_base_idx,
                           esp.nv,
                           num_cols,
                           "layer",
                           "delta_col_mass",
                           alf.delta_col_mass,
                           "deltacolmassprofile",
                           false);
    }
    else {
        print_data_to_file(esp,
                           nstep,
                           column_base_idx,
                           esp.nv,
                           num_cols,
                           "layer",
                           "delta_col_upper",
                           alf.delta_col_upper,
                           "deltacolupperprofile",
                           false);

        print_data_to_file(esp,
                           nstep,
                           column_base_idx,
                           esp.nv,
                           num_cols,
                           "layer",
                           "delta_col_lower",
                           alf.delta_col_lower,
                           "deltacollowerprofile",
                           false);
    }

    {
        // Print out opacities data
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.opac_wg_lay,
                                         "opacprofile",
                                         false);
    }

    if (!iso) {
        // Print out opacities data
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nvi,
                                         num_cols,
                                         "interface",
                                         alf.opac_wg_int,
                                         "opacintprofile",
                                         false);
    }

    //***********************************************************************************************
    if (iso) {
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.delta_tau_wg,
                                         "opt_depthprofile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.trans_wg,
                                         "trans_band_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.w0_wg,
                                         "single_scat_band_profile",
                                         false);
        print_weighted_band_data_to_file(
            esp, nstep, column_base_idx, esp.nv, num_cols, "layer", alf.M_term, "M_profile", false);
        print_weighted_band_data_to_file(
            esp, nstep, column_base_idx, esp.nv, num_cols, "layer", alf.N_term, "N_profile", false);
        print_weighted_band_data_to_file(
            esp, nstep, column_base_idx, esp.nv, num_cols, "layer", alf.P_term, "P_profile", false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_plus,

                                         "G_plus_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_minus,
                                         "G_minus_profile",
                                         false);
    }
    else {
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.delta_tau_wg_upper,

                                         "opt_depth_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.delta_tau_wg_lower,

                                         "opt_depth_lower_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.trans_wg_upper,
                                         "trans_band_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.trans_wg_lower,
                                         "trans_band_lower_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.w0_wg_upper,
                                         "single_scat_band_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.w0_wg_lower,
                                         "single_scat_band_lower_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.M_upper,
                                         "M_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.M_lower,
                                         "M_lower_profile",
                                         false);

        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.N_upper,
                                         "N_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.N_lower,
                                         "N_lower_profile",
                                         false);

        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.P_upper,
                                         "P_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.P_lower,
                                         "P_lower_profile",
                                         false);

        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_plus_upper,
                                         "G_plus_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_plus_lower,
                                         "G_plus_lower_profile",
                                         false);

        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_minus_upper,
                                         "G_minus_upper_profile",
                                         false);
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nv,
                                         num_cols,
                                         "layer",
                                         alf.G_minus_lower,
                                         "G_minus_lower_profile",
                                         false);
    }
    //***********************************************************************************************
    {
        // Print out downward flux
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nvi,
                                         num_cols,
                                         "interface",
                                         F_down_wg,
                                         "F_down_profile",
                                         false);
    }

    {
        // Print out downward flux
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nvi,
                                         num_cols,
                                         "interface",
                                         F_up_wg,
                                         "F_up_profile",
                                         false);
    }

    {
        // Print out direct beam flux
        print_weighted_band_data_to_file(esp,
                                         nstep,
                                         column_base_idx,
                                         esp.nvi,
                                         num_cols,
                                         "interface",
                                         F_dir_wg,
                                         "F_dir_profile",
                                         false);
    }

    if (!iso) {
        if (thomas) {
            // Print out thomas sol (gives interface sol
            print_X_buff_thomas_data_to_file(esp,
                                             nstep,
                                             column_base_idx,
                                             2 * (esp.nv * 2 + 1),
                                             num_cols,
                                             "int",
                                             alf.X_buff,
                                             "X_thomas_profile");
        }

        {
            // Print out direct beam flux
            print_weighted_band_data_to_file(esp,
                                             nstep,
                                             column_base_idx,
                                             esp.nv,
                                             num_cols,
                                             "layer",
                                             Fc_dir_wg,
                                             "Fc_dir_profile",
                                             false);
        }
    }

    {
        // Print out alf qheat

        for (int c = 0; c < num_cols; c++) {
            int col_offset = (column_base_idx + c) * esp.nv;

            std::string DBG_OUTPUT_DIR = esp.get_output_dir()
                                         + "/alfprof/"
                                           "step_"
                                         + std::to_string(nstep) + "/column_"
                                         + std::to_string(column_base_idx + c) + "/";
            create_output_dir(DBG_OUTPUT_DIR);

            string output_file_name = DBG_OUTPUT_DIR + "alf_qheat_profile.dat";

            FILE *output_file = fopen(output_file_name.c_str(), "w");

            std::shared_ptr<double[]> qh_h = get_cuda_data(&((Qheat.ptr())[col_offset]), esp.nv);

            fprintf(output_file, "level\tqheat\n");

            for (int i = 0; i < esp.nv; i++) {
                fprintf(output_file, "%d\t%#.6g\n", i, qh_h[i]);
            }

            fclose(output_file);
            cuda_check_status_or_exit(string(__FILE__ ":"
                                                      "qheat")
                                          .c_str(),
                                      __LINE__);
        }
    }
}
