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

        if (coszrs < 0.0)
            cos_zenith_angles[column_idx] = 0.0;
        else
            cos_zenith_angles[column_idx] = coszrs;
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

    return true;
}

bool Insolation::initialise_memory(const ESP&               esp,
                                   device_RK_array_manager& phy_modules_core_arrays) {
    if (enabled) {
        cos_zenith_angles.allocate(esp.point_num);
        cos_zenith_angles.zero();
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
}

bool Insolation::initial_conditions(const ESP& esp, const SimulationSetup& sim, storage* s) {
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
    }

    return true;
}

bool Insolation::phy_loop(ESP&                   esp,
                          const SimulationSetup& sim,
                          int                    nstep, // Step number
                          double                 time_step) {

    const int num_blocks = 256;
    if (enabled) {
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
