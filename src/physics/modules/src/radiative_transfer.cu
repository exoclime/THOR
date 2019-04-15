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
// ESP -  Exoclimes Simulation Platform. (version 1.0)
//
//
//
// Method: Radiative transfer physics module
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#include "radiative_transfer.h"

#include "profx_RT.h"

radiative_transfer::radiative_transfer() {
}

radiative_transfer::~radiative_transfer() {
}

void radiative_transfer::print_config() {
    log::printf("  Radiative transfer module\n");

    // basic star-planet properties
    log::printf("    Tstar                       = %f K.\n", Tstar_config);
    log::printf("    Orbital distance            = %f au.\n", planet_star_dist_config);
    log::printf("    Radius of host star         = %f R_sun.\n", radius_star_config);
    log::printf("    Diffusivity factor          = %f.\n", diff_fac_config);
    log::printf("    Lower boundary temperature  = %f K.\n", Tlow_config);
    log::printf("    Bond albedo                 = %f.\n", albedo_config);
    log::printf("    Shortwave Absorption coef   = %f.\n", tausw_config);
    log::printf("    Longwave Absorption coef    = %f.\n", taulw_config);
    log::printf("    Using sin(lat) variation LW?    = %f.\n", latf_lw_config);
    log::printf("    Longwave Absorption coef (poles)    = %f.\n", taulw_pole_config);
    log::printf("\n");

    // orbit/insolation properties
    log::printf("    Synchronous rotation        = %s.\n", sync_rot_config ? "true" : "false");
    log::printf("    Orbital mean motion         = %f rad/s.\n", mean_motion_config);
    log::printf("    Host star initial right asc = %f deg.\n", alpha_i_config);
    log::printf("    Planet true initial long    = %f.\n", true_long_i_config);
    log::printf("    Orbital eccentricity        = %f.\n", ecc_config);
    log::printf("    Obliquity                   = %f deg.\n", obliquity_config);
    log::printf("    Longitude of periastron     = %f deg.\n", longp_config);

    // surface parameters
    log::printf("    Surface                     = %s.\n", surface_config ? "true" : "false");
    log::printf("    Surface Heat Capacity       = %f deg.\n", Csurf_config);
}

bool radiative_transfer::initialise_memory(const ESP &              esp,
                                           device_RK_array_manager &phy_modules_core_arrays) {
    //  Rad Transfer
    cudaMalloc((void **)&fnet_up_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&fnet_dn_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&tau_d, esp.nv * esp.point_num * 2 * sizeof(double));

    cudaMalloc((void **)&phtemp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&thtemp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&ttemp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dtemp, esp.nv * esp.point_num * sizeof(double));

    insol_h = (double *)malloc(esp.point_num * sizeof(double));
    cudaMalloc((void **)&insol_d, esp.point_num * sizeof(double));
    insol_ann_h = (double *)malloc(esp.point_num * sizeof(double));
    cudaMalloc((void **)&insol_ann_d, esp.point_num * sizeof(double));

    fnet_up_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    fnet_dn_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    tau_h     = (double *)malloc(esp.nv * esp.point_num * 2 * sizeof(double));

    cudaMalloc((void **)&surf_flux_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&Tsurface_d, esp.point_num * sizeof(double));

    return true;
}


bool radiative_transfer::free_memory() {
    cudaFree(fnet_up_d);
    cudaFree(fnet_dn_d);
    cudaFree(tau_d);

    cudaFree(phtemp);
    cudaFree(thtemp);
    cudaFree(ttemp);
    cudaFree(dtemp);

    return true;
}

bool radiative_transfer::initial_conditions(const ESP &esp, const SimulationSetup &sim) {
    RTSetup(Tstar_config,
            planet_star_dist_config,
            radius_star_config,
            diff_fac_config,
            Tlow_config,
            albedo_config,
            tausw_config,
            taulw_config,
            latf_lw_config,
            taulw_pole_config,
            sync_rot_config,
            mean_motion_config,
            true_long_i_config,
            longp_config,
            ecc_config,
            alpha_i_config,
            obliquity_config,
            sim.Omega,
            surface_config,
            Csurf_config,
            esp.point_num);
    return true;
}

bool radiative_transfer::phy_loop(ESP &                  esp,
                                  const SimulationSetup &sim,
                                  int                    nstep, // Step number
                                  double                 time_step) {           // Time-step [s]

    //  update global insolation properties if necessary
    if (sync_rot) {
        if (ecc > 1e-10) {
            update_spin_orbit(nstep * time_step, sim.Omega);
        }
    }
    else {
        update_spin_orbit(nstep * time_step, sim.Omega);
    }

    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3 NB((esp.point_num / NTH) + 1, esp.nv, 1);
    dim3 NBRT((esp.point_num / NTH) + 1, 1, 1);

    rtm_dual_band<<<NBRT, NTH>>>(esp.pressure_d,
                                 //rtm_dual_band <<< 1,1 >>> (pressure_d         ,
                                 esp.Rho_d,
                                 esp.temperature_d,
                                 fnet_up_d,
                                 fnet_dn_d,
                                 tau_d,
                                 sim.Gravit,
                                 sim.Cp,
                                 esp.lonlat_d,
                                 esp.Altitude_d,
                                 esp.Altitudeh_d,
                                 phtemp,
                                 dtemp,
                                 ttemp,
                                 thtemp,
                                 time_step,
                                 Tstar,
                                 planet_star_dist,
                                 radius_star,
                                 diff_fac,
                                 Tlow,
                                 albedo,
                                 tausw,
                                 taulw,
                                 latf_lw,
                                 taulw_pole,
                                 incflx,
                                 sim.P_Ref,
                                 esp.point_num,
                                 esp.nv,
                                 esp.nvi,
                                 sim.A,
                                 r_orb,
                                 alpha, //current RA of star (relative to zero long on planet)
                                 alpha_i,
                                 sin_decl, //declination of star
                                 cos_decl,
                                 sync_rot,
                                 ecc,
                                 obliquity,
                                 insol_d,
                                 surface,
                                 Csurf,
                                 Tsurface_d,
                                 surf_flux_d);

    if (nstep * time_step < (2 * M_PI / mean_motion)) {
        // stationary orbit/obliquity
        // calculate annually average of insolation for the first orbit
        annual_insol<<<NBRT, NTH>>>(insol_ann_d, insol_d, nstep);
    }
    return true;
}

bool radiative_transfer::configure(config_file &config_reader) {
    // basic star-planet properties
    config_reader.append_config_var("Tstar", Tstar_config, Tstar_config);
    config_reader.append_config_var(
        "planet_star_dist", planet_star_dist_config, planet_star_dist_config);
    config_reader.append_config_var("radius_star", radius_star_config, radius_star_config);
    config_reader.append_config_var("diff_fac", diff_fac_config, diff_fac_config);
    config_reader.append_config_var("Tlow", Tlow_config, Tlow_config);
    config_reader.append_config_var("albedo", albedo_config, albedo_config);
    config_reader.append_config_var("tausw", tausw_config, tausw_config);
    config_reader.append_config_var("taulw", taulw_config, taulw_config);

    // options for latitude dependence in longwave opacity
    config_reader.append_config_var("latf_lw", latf_lw_config, latf_lw_config);
    config_reader.append_config_var("taulw_pole", taulw_pole_config, taulw_pole_config);

    // orbit/insolation properties
    config_reader.append_config_var("sync_rot", sync_rot_config, sync_rot_config);
    config_reader.append_config_var("mean_motion", mean_motion_config, mean_motion_config);
    config_reader.append_config_var("alpha_i", alpha_i_config, alpha_i_config);
    config_reader.append_config_var("true_long_i", true_long_i_config, true_long_i_config);
    config_reader.append_config_var("ecc", ecc_config, ecc_config);
    config_reader.append_config_var("obliquity", obliquity_config, obliquity_config);
    config_reader.append_config_var("longp", longp_config, longp_config);

    // properties for a solid/liquid surface
    config_reader.append_config_var("surface", surface_config, surface_config);
    config_reader.append_config_var("Csurf", Csurf_config, Csurf_config);


    return true;
}

bool radiative_transfer::store(const ESP &esp, storage &s) {
    cudaMemcpy(insol_h, insol_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(insol_h, esp.point_num, "/insol", "W m^-2", "insolation (instantaneous)");

    cudaMemcpy(insol_ann_h, insol_ann_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(insol_ann_h,
                   esp.point_num,
                   "/insol_annual",
                   "W m^-2",
                   "insolation (annual/orbit averaged)");

    cudaMemcpy(
        fnet_up_h, fnet_up_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(fnet_up_h, esp.nvi * esp.point_num, "/fnet_up", "W m^-2", "upward flux");

    cudaMemcpy(
        fnet_dn_h, fnet_dn_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(fnet_dn_h, esp.nvi * esp.point_num, "/fnet_dn", "W m^-2", "downward flux");

    cudaMemcpy(tau_h, tau_d, esp.nv * esp.point_num * 2 * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(
        tau_h, esp.nv * esp.point_num * 2, "/tau", " ", "optical depth across each layer");

    return true;
}

bool radiative_transfer::store_init(storage &s) {
    s.append_value(Tstar, "/Tstar", "K", "Temperature of host star");
    s.append_value(Tlow, "/Tlow", "K", "Temperature of interior heat flux");
    s.append_value(planet_star_dist / 149597870.7,
                   "/planet_star_dist",
                   "au",
                   "distance b/w host star and planet");
    s.append_value(radius_star / 695508, "/radius_star", "R_sun", "radius of host star");
    s.append_value(diff_fac, "/diff_fac", "-", "diffusivity factor");
    s.append_value(albedo, "/albedo", "-", "bond albedo of planet");
    s.append_value(tausw, "/tausw", "-", "shortwave optical depth of deepest layer");
    s.append_value(taulw, "/taulw", "-", "longwave optical depth of deepest layer");
    s.append_value(sync_rot ? 1.0 : 0.0, "/sync_rot", "-", "enforce synchronous rotation");
    s.append_value(mean_motion, "/mean_motion", "rad/s", "Orbital mean motion");
    s.append_value(alpha_i * 180 / M_PI, "/alpha_i", "deg", "initial RA of host star");
    s.append_value(
        true_long_i * 180 / M_PI, "/true_long_i", "deg", "initial orbital position of planet");
    s.append_value(ecc, "/ecc", "-", "orbital eccentricity");
    s.append_value(obliquity * 180 / M_PI, "/obliquity", "deg", "tilt of spin axis");
    s.append_value(longp * 180 / M_PI, "/longp", "deg", "longitude of periastron");

    s.append_value(latf_lw ? 1.0 : 0.0, "/latf_lw", "-", "use lat dependent opacity");
    s.append_value(
        taulw_pole, "/taulw_pole", "-", "longwave optical depth of deepest layer at poles");

    s.append_value(surface ? 1.0 : 0.0, "/surface", "-", "include solid/liquid surface");
    s.append_value(Csurf, "/Csurf", "J/K/m^2", "heat capacity of surface by area");


    return true;
}

void radiative_transfer::RTSetup(double Tstar_,
                                 double planet_star_dist_,
                                 double radius_star_,
                                 double diff_fac_,
                                 double Tlow_,
                                 double albedo_,
                                 double tausw_,
                                 double taulw_,
                                 bool   latf_lw_,
                                 double taulw_pole_,
                                 bool   sync_rot_,
                                 double mean_motion_,
                                 double true_long_i_,
                                 double longp_,
                                 double ecc_,
                                 double alpha_i_,
                                 double obliquity_,
                                 double Omega,
                                 bool   surface_,
                                 double Csurf_,
                                 int    point_num) {

    double bc = 5.677036E-8; // Stefan–Boltzmann constant [W m−2 K−4]

    Tstar            = Tstar_;
    planet_star_dist = planet_star_dist_ * 149597870.7; //conv to km
    radius_star      = radius_star_ * 695508;           //conv to km
    diff_fac         = diff_fac_;
    Tlow             = Tlow_;
    albedo           = albedo_;
    tausw            = tausw_;
    taulw            = taulw_;
    taulw_pole       = taulw_pole_;
    latf_lw          = latf_lw_;
    double resc_flx  = pow(radius_star / planet_star_dist, 2.0);
    incflx           = resc_flx * bc * Tstar * Tstar * Tstar * Tstar;

    sync_rot = sync_rot_;
    if (sync_rot) {
        mean_motion = Omega; //just set for the sync_rot, obl != 0 case
    }
    else {
        mean_motion = mean_motion_;
    }
    true_long_i           = true_long_i_ * M_PI / 180.0;
    longp                 = longp_ * M_PI / 180.0;
    ecc                   = ecc_;
    double true_anomaly_i = fmod(true_long_i - longp, (2 * M_PI));
    double ecc_anomaly_i  = true2ecc_anomaly(true_anomaly_i, ecc);
    mean_anomaly_i        = fmod(ecc_anomaly_i - ecc * sin(ecc_anomaly_i), (2 * M_PI));
    alpha_i               = alpha_i_ * M_PI / 180.0;
    obliquity             = obliquity_ * M_PI / 180.0;

    surface = surface_;
    Csurf   = Csurf_;
}

void radiative_transfer::update_spin_orbit(double time, double Omega) {

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
