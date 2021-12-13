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
//
// Description: Insolation computation
//
//
//
// Known limitations: None
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
#include "insolation.h"

#include "debug_helpers.h"

#include "binary_test.h"
#include "debug.h"


// ***************************************************************************************************
//** Kernels doing the real computation work

__host__ double sign(double val) {
    if (val < 0.0)
        return -1.0;
    else
        return 1.0;
}

__host__ double solve_kepler(double mean_anomaly, double ecc) {

    // Solve Kepler's equation (see Murray & Dermott 1999)
    // Get eccentric anomaly from mean anomaly and eccentricity

    double ecc_anomaly, fi, fi_1, fi_2, fi_3, di_1, di_2, di_3;

    ecc_anomaly = mean_anomaly + sign(sin(mean_anomaly)) * 0.85 * ecc;
    di_3        = 1.0;

    while (di_3 > 1e-15) {
        fi          = ecc_anomaly - ecc * sin(ecc_anomaly) - mean_anomaly;
        fi_1        = 1.0 - ecc * cos(ecc_anomaly);
        fi_2        = ecc * sin(ecc_anomaly);
        fi_3        = ecc * cos(ecc_anomaly);
        di_1        = -fi / fi_1;
        di_2        = -fi / (fi_1 + 0.5 * di_1 * fi_2);
        di_3        = -fi / (fi_1 + 0.5 * di_2 * fi_2 + 1. / 6. * di_2 * di_2 * fi_3);
        ecc_anomaly = ecc_anomaly + di_3;
    }
    return ecc_anomaly;
}

__host__ double calc_r_orb(double ecc_anomaly, double ecc) {

    // Calc relative distance between planet and star (units of semi-major axis)

    double r = 1.0 - ecc * cos(ecc_anomaly);
    return r;
}

__host__ double ecc2true_anomaly(double ecc_anomaly, double ecc) {

    // Convert from eccentric to true anomaly

    double tanf2, true_anomaly;
    tanf2        = sqrt((1.0 + ecc) / (1.0 - ecc)) * tan(ecc_anomaly / 2.0);
    true_anomaly = 2.0 * atan(tanf2);
    if (true_anomaly < 0.0)
        true_anomaly += 2 * M_PI;
    return true_anomaly;
}

__host__ double true2ecc_anomaly(double true_anomaly, double ecc) {

    // Convert from true to eccentric anomaly

    double cosE, ecc_anomaly;
    while (true_anomaly < 0.0)
        true_anomaly += 2 * M_PI;
    while (true_anomaly >= 2 * M_PI)
        true_anomaly -= 2 * M_PI;
    cosE = (cos(true_anomaly) + ecc) / (1.0 + ecc * cos(true_anomaly));
    if (true_anomaly < M_PI) {
        ecc_anomaly = acos(cosE);
    }
    else {
        ecc_anomaly = 2 * M_PI - acos(cosE);
    }
    return ecc_anomaly;
}

__device__ double
calc_zenith(double*      lonlat_d, //latitude/longitude grid
            const double alpha,    //current RA of star (relative to zero long on planet)
            const double alpha_i,
            const double sin_decl, //declination of star
            const double cos_decl,
            const bool   sync_rot,
            const double ecc,
            const double obliquity,
            const int    id) {

    // Calculate the insolation (scaling) at a point on the surface

    double coszrs;

    if (sync_rot) {
        if (ecc < 1e-10) {
            if (obliquity < 1e-10) { //standard sync, circular, zero obl case
                coszrs = cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2] - alpha_i);
            }
            else { //sync, circular, but some obliquity
                coszrs = (sin(lonlat_d[id * 2 + 1]) * sin_decl
                          + cos(lonlat_d[id * 2 + 1]) * cos_decl * cos(lonlat_d[id * 2] - alpha_i));
            }
        }
        else {                       //in below cases, watch out for numerical drift of mean(alpha)
            if (obliquity < 1e-10) { // sync, zero obliquity, but ecc orbit
                coszrs = cos(lonlat_d[id * 2 + 1]) * cos(lonlat_d[id * 2] - alpha);
            }
            else { // sync, non-zero obliquity, ecc orbit (full calculation applies)
                coszrs = (sin(lonlat_d[id * 2 + 1]) * sin_decl
                          + cos(lonlat_d[id * 2 + 1]) * cos_decl * cos(lonlat_d[id * 2] - alpha));
            }
        }
    }
    else {
        coszrs = (sin(lonlat_d[id * 2 + 1]) * sin_decl
                  + cos(lonlat_d[id * 2 + 1]) * cos_decl * cos(lonlat_d[id * 2] - alpha));
    }
    return coszrs; //zenith angle
}

// compute zenith angle for a full grid of lat/lon data
__global__ void compute_cos_zenith_angles(double* cos_zenith_angles,
                                          double* lonlat_d,
                                          double  alpha,
                                          double  alpha_i,
                                          double  sin_decl,
                                          double  cos_decl,
                                          double  ecc,
                                          double  obliquity,
                                          bool    sync_rot,
                                          int     num_points) {
    // helios_angle_star = pi - zenith_angle
    // cos(helios_angle_star) = mu_star = cos(pi - zenith_angle) = -cos(zenith_angle)
    // Zenith angle is only positive.
    // mu_star is only negative in helios -> need to change sign outside of here for Alfrodull
    // radiative transfer uses zenith angle
    int column_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (column_idx < num_points) {
        double coszrs = calc_zenith(lonlat_d, //latitude/longitude grid
                                    alpha,    //current RA of star (relative to zero long on planet)
                                    alpha_i,
                                    sin_decl, //declination of star
                                    cos_decl,
                                    sync_rot,
                                    ecc,
                                    obliquity,
                                    column_idx);

        //hack
        // double coszrs = cos(0 * M_PI / 180.);

        if (coszrs < 0.0)
            cos_zenith_angles[column_idx] = 0.0;
        else
            cos_zenith_angles[column_idx] = coszrs;
    }
}

// add to zenith angle average (trapz rule)
__global__ void trapz_step_daily_avg(double* cos_zenith_daily,
                                     double* cos_zenith_angles,
                                     double  dtstep,
                                     double  Pday,
                                     int     iday,
                                     int     istep,
                                     int     n_day_steps,
                                     int     num_points) {

    int column_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (column_idx < num_points) {
        if (istep == 0 || istep == n_day_steps) { //endpoints of trapz integral
            cos_zenith_daily[iday * num_points + column_idx] +=
                0.5 * dtstep * cos_zenith_angles[column_idx] / Pday;
        }
        else { //interior points
            cos_zenith_daily[iday * num_points + column_idx] +=
                dtstep * cos_zenith_angles[column_idx] / Pday;
        }
    }
}


// ***************************************************************************************************
// ** Insolation class doing the management stuff
Insolation::Insolation() {
}

bool Insolation::configure(config_file& config_reader) {
    // orbit/insolation properties
    config_reader.append_config_var("sync_rot", sync_rot_config, sync_rot_config);
    config_reader.append_config_var("mean_motion", mean_motion_config, mean_motion_config);
    config_reader.append_config_var("alpha_i", alpha_i_config, alpha_i_config);
    config_reader.append_config_var("true_long_i", true_long_i_config, true_long_i_config);
    config_reader.append_config_var("ecc", ecc_config, ecc_config);
    config_reader.append_config_var("obliquity", obliquity_config, obliquity_config);
    config_reader.append_config_var("longp", longp_config, longp_config);

    config_reader.append_config_var("insol_avg", insol_avg_str, string(insol_avg_default));
    return true;
}

bool Insolation::initialise_memory(const ESP&               esp,
                                   device_RK_array_manager& phy_modules_core_arrays) {
    if (enabled) {
        cos_zenith_angles.allocate(esp.point_num);
        cos_zenith_angles.zero();

        USE_BENCHMARK();

#ifdef BENCHMARKING
        std::map<string, output_def> debug_arrays = {
            {"coszs",
             {cos_zenith_angles.ptr_ref(),
              (int)cos_zenith_angles.get_size(),
              "coszs",
              "cz",
              true,
              dummy}},
        };

        BENCH_POINT_REGISTER_PHY_VARS(debug_arrays, (), ());
#endif // BENCHMARKING
    }

    return true;
}

void Insolation::print_config() {
    log::printf("  Insolation configuration (if used by modules) \n");

    // orbit/insolation properties
    log::printf("    Synchronous rotation        = %s.\n", sync_rot_config ? "true" : "false");
    log::printf("    Orbital mean motion         = %f rad/s.\n", mean_motion_config);
    log::printf("    Host star initial right asc = %f deg.\n", alpha_i_config);
    log::printf("    Planet true initial long    = %f.\n", true_long_i_config);
    log::printf("    Orbital eccentricity        = %f.\n", ecc_config);
    log::printf("    Obliquity                   = %f deg.\n", obliquity_config);
    log::printf("    Longitude of periastron     = %f deg.\n", longp_config);
    log::printf("    Use averaged insolation     = %s \n", insol_avg_str.c_str());
}

bool Insolation::initial_conditions(const ESP& esp, const SimulationSetup& sim, storage* s) {
    bool config_OK = true;
    if (enabled) {
        sync_rot = sync_rot_config;
        if (sync_rot) {
            mean_motion = sim.Omega; //just set for the sync_rot, obl != 0 case
        }
        else {
            mean_motion = mean_motion_config;
        }

        true_long_i           = true_long_i_config * M_PI / 180.0;
        longp                 = longp_config * M_PI / 180.0;
        ecc                   = ecc_config;
        double true_anomaly_i = fmod(true_long_i - longp, (2 * M_PI));
        double ecc_anomaly_i  = true2ecc_anomaly(true_anomaly_i, ecc);
        mean_anomaly_i        = fmod(ecc_anomaly_i - ecc * sin(ecc_anomaly_i), (2 * M_PI));
        alpha_i               = alpha_i_config * M_PI / 180.0;
        obliquity             = obliquity_config * M_PI / 180.0;

        insol_avg = NO_INSOL_AVG;
        if (insol_avg_str == "NoInsolAvg") {
            insol_avg = NO_INSOL_AVG;
            config_OK &= true;
        }
        else if (insol_avg_str == "DiurnalAvg" || insol_avg_str == "Diurnal") {
            //average over mean "day" (Pday = 2 PI/(Omega - n_orb))
            if (sync_rot) {
                printf("Diurnal average insolation cannot be used with synchronous rotation\n");
                config_OK &= false;
            }
            else {
                insol_avg = DIURNAL_AVG;
                config_OK &= true;
            }
        }
        else if (insol_avg_str == "AnnualAvg" || insol_avg_str == "Annual") {
            //average over orbit (Porb = 2 PI/n_orb)
            insol_avg = ANNUAL_AVG;
            printf("Annual average insolation option not ready yet...\n");
            config_OK &= false; //not ready yet!!
        }
        else {
            log::printf("insol_avg config item not recognised: [%s]\n", insol_avg_str.c_str());
            config_OK &= false;
        }

        // compute averages if necessary---------------------------------------
        const int num_blocks = 256;
        if (insol_avg == DIURNAL_AVG) {
            if (sim.Omega != mean_motion) {
                Pday = 2 * M_PI / (sim.Omega - mean_motion); //what do i do in case of Omega < 0 ??
            }
            else {
                printf("Rotation and orbital periods cannot be equal with diurnal avg forcing\n");
                config_OK &= false;
            }
            if (config_OK) {
                Porb       = 2 * M_PI / mean_motion;
                n_days_orb = int(Porb / Pday); //hmm what to do with remaining fraction of day??
                cos_zenith_daily.allocate(n_days_orb * esp.point_num);
                cos_zenith_daily.zero();
                day_start_time.allocate(n_days_orb);
                day_start_time.zero();

                std::shared_ptr<double[]> day_start_time_h = day_start_time.get_host_data_ptr();
                double                    dtstep           = Pday / n_day_steps;
                for (int iday = 0; iday < n_days_orb; iday++) {
                    //loop over days and compute avg insolation/zenith angle
                    //compute day start time
                    day_start_time_h[iday] = iday * Pday;

                    for (int istep = 0; istep <= n_day_steps; istep++) {
                        //trapezoidal rule here
                        update_spin_orbit(day_start_time_h[iday] + istep * dtstep, sim.Omega);

                        compute_cos_zenith_angles<<<(esp.point_num / num_blocks) + 1, num_blocks>>>(
                            *cos_zenith_angles,
                            esp.lonlat_d,
                            alpha,
                            alpha_i,
                            sin_decl,
                            cos_decl,
                            ecc,
                            obliquity,
                            sync_rot,
                            esp.point_num);
                        cudaDeviceSynchronize();
                        cuda_check_status_or_exit(__FILE__, __LINE__);

                        trapz_step_daily_avg<<<(esp.point_num / num_blocks) + 1, num_blocks>>>(
                            *cos_zenith_daily,
                            *cos_zenith_angles,
                            dtstep,
                            Pday,
                            iday,
                            istep,
                            n_day_steps,
                            esp.point_num);

                        cudaDeviceSynchronize();
                        cuda_check_status_or_exit(__FILE__, __LINE__);
                    }
                    // printf("day = %d, time = %f, M = %f\n",
                    //        iday,
                    //        day_start_time_h[iday],
                    //        mean_anomaly);
                }
                //reset orbit to original position
                update_spin_orbit(0.0, sim.Omega);

                day_start_time.put();
            }
        }

        if (!config_OK) {
            log::printf("Error in configuration file\n");
            exit(-1);
        }
    }

    return true;
}

bool Insolation::phy_loop(ESP&                   esp,
                          const SimulationSetup& sim,
                          int                    nstep, // Step number
                          double                 time_step) {

    const int num_blocks = 256;
    if (enabled) {
        USE_BENCHMARK();

        //  update global insolation properties if necessary
        if (sync_rot) {
            if (ecc > 1e-10) {
                update_spin_orbit(nstep * time_step, sim.Omega);
            }
        }
        else {
            update_spin_orbit(nstep * time_step, sim.Omega);
        }


        compute_cos_zenith_angles<<<(esp.point_num / num_blocks) + 1, num_blocks>>>(
            *cos_zenith_angles,
            esp.lonlat_d,
            alpha,
            alpha_i,
            sin_decl,
            cos_decl,
            ecc,
            obliquity,
            sync_rot,
            esp.point_num);
        cudaDeviceSynchronize();
        cuda_check_status_or_exit(__FILE__, __LINE__);

        BENCH_POINT_I_PHY(nstep, "Insolation", (), ("coszs"));
    }


    return true;
}


bool Insolation::store_init(storage& s) {
    if (enabled) {
        if (!s.has_table("/sync_rot"))
            s.append_value(sync_rot ? 1.0 : 0.0, "/sync_rot", "-", "enforce synchronous rotation");
        if (!s.has_table("/mean_motion"))
            s.append_value(mean_motion, "/mean_motion", "rad/s", "orbital mean motion");
        if (!s.has_table("/alpha_i"))
            s.append_value(alpha_i * 180 / M_PI, "/alpha_i", "deg", "initial RA of host star");
        if (!s.has_table("/true_long_i"))
            s.append_value(true_long_i * 180 / M_PI,
                           "/true_long_i",
                           "deg",
                           "initial orbital position of planet");
        if (!s.has_table("/ecc"))
            s.append_value(ecc, "/ecc", "-", "orbital eccentricity");
        if (!s.has_table("/obliquity"))
            s.append_value(obliquity * 180 / M_PI, "/obliquity", "deg", "tilt of spin axis");
        if (!s.has_table("/longp"))
            s.append_value(longp * 180 / M_PI, "/longp", "deg", "longitude of periastron");
    }


    return true;
}

bool Insolation::store(const ESP& esp, storage& s) {
    if (enabled) {
        std::shared_ptr<double[]> cos_zenith_angles_h = cos_zenith_angles.get_host_data();
        s.append_table(cos_zenith_angles_h.get(),
                       cos_zenith_angles.get_size(),
                       "/cos_zenith_angles",
                       "-",
                       "cos_zenith_angles per column");

        if (insol_avg == DIURNAL_AVG) {
            std::shared_ptr<double[]> cos_zenith_daily_h = cos_zenith_daily.get_host_data();
            s.append_table(cos_zenith_daily_h.get(),
                           cos_zenith_daily.get_size(),
                           "/cos_zenith_daily",
                           "-",
                           "cos_zenith (daily avg) per column");

            std::shared_ptr<double[]> day_start_time_h = day_start_time.get_host_data();
            s.append_table(day_start_time_h.get(),
                           day_start_time.get_size(),
                           "/day_start_time",
                           "s",
                           "start time of each day");
        }
    }

    return true;
}


void Insolation::update_spin_orbit(double time, double Omega) {

    // Update the insolation related parameters for spin and orbit
    double ecc_anomaly, true_long;

    mean_anomaly = fmod((mean_motion * time + mean_anomaly_i), (2 * M_PI));

    ecc_anomaly = fmod(solve_kepler(mean_anomaly, ecc), (2 * M_PI));

    r_orb     = calc_r_orb(ecc_anomaly, ecc);
    true_long = fmod((ecc2true_anomaly(ecc_anomaly, ecc) + longp), (2 * M_PI));

    sin_decl = sin(obliquity) * sin(true_long);
    cos_decl = sqrt(1.0 - sin_decl * sin_decl);
    alpha    = -Omega * time + true_long - true_long_i + alpha_i;
}
