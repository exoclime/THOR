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

#include "insolation.h"

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
    log::printf("    1.0/Diffusivity factor      = %f.\n", diff_ang_config);
    // log::printf("    Internal flux temperature   = %f K.\n", Tint_config);
    log::printf("    Bond albedo                 = %f.\n", albedo_config);
    // log::printf("    Shortwave Absorption coef   = %f.\n", tausw_config);
    // log::printf("    Longwave Absorption coef    = %f.\n", taulw_config);
    log::printf("    Using sin(lat) variation LW?      = %s.\n", latf_lw_config ? "true" : "false");
    log::printf("    Longwave opacity (poles)  = %f.\n", kappa_lw_pole_config);
    log::printf("    Power law index of unmixed LW abs = %f.\n", n_lw_config);
    // log::printf("    Strength of mixed LW abs    = %f.\n", f_lw_config);
    log::printf("    Power law index of SW       = %f.\n", n_sw_config);
    log::printf("\n");


    log::printf("    1D mode                     = %s.\n", rt1Dmode_config ? "true" : "false");

    // spinup-spindown parameters
    log::printf("    Spin up start step          = %d.\n", spinup_start_step);
    log::printf("    Spin up stop step           = %d.\n", spinup_stop_step);
    log::printf("    Spin down start step        = %d.\n", spindown_start_step);
    log::printf("    Spin down stop step         = %d.\n", spindown_stop_step);
}

bool radiative_transfer::initialise_memory(const ESP &              esp,
                                           device_RK_array_manager &phy_modules_core_arrays) {
    //  Rad Transfer
    cudaMalloc((void **)&flw_up_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&flw_dn_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&fsw_up_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&fsw_dn_d, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&tau_d, esp.nv * esp.point_num * 2 * sizeof(double));

    cudaMalloc((void **)&phtemp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&thtemp, esp.nvi * esp.point_num * sizeof(double));
    cudaMalloc((void **)&ttemp, esp.nv * esp.point_num * sizeof(double));
    cudaMalloc((void **)&dtemp, esp.nv * esp.point_num * sizeof(double));

    cudaMalloc((void **)&qheat_d, esp.nv * esp.point_num * sizeof(double));
    qheat_h = (double *)malloc(esp.point_num * esp.nv * sizeof(double));

    insol_h = (double *)malloc(esp.point_num * sizeof(double));
    cudaMalloc((void **)&insol_d, esp.point_num * sizeof(double));
    insol_ann_h = (double *)malloc(esp.point_num * sizeof(double));
    cudaMalloc((void **)&insol_ann_d, esp.point_num * sizeof(double));

    flw_up_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    flw_dn_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    tau_h    = (double *)malloc(esp.nv * esp.point_num * 2 * sizeof(double));

    fsw_up_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
    fsw_dn_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));

    cudaMalloc((void **)&surf_flux_d, esp.point_num * sizeof(double));

    return true;
}


bool radiative_transfer::free_memory() {
    cudaFree(flw_up_d);
    cudaFree(flw_dn_d);
    cudaFree(fsw_up_d);
    cudaFree(fsw_dn_d);
    cudaFree(tau_d);

    cudaFree(phtemp);
    cudaFree(thtemp);
    cudaFree(ttemp);
    cudaFree(dtemp);

    cudaFree(qheat_d);
    free(qheat_h);

    cudaFree(insol_d);
    cudaFree(insol_ann_d);

    free(flw_up_h);
    free(flw_dn_h);
    free(tau_h);

    free(fsw_up_h);
    free(fsw_dn_h);

    cudaFree(surf_flux_d);

    return true;
}

bool radiative_transfer::initial_conditions(const ESP &            esp,
                                            const SimulationSetup &sim,
                                            storage *              s) {

    if (spinup_start_step > -1 || spinup_stop_step > -1) {
        if (spinup_stop_step < spinup_start_step)
            printf("DGRT: inconsistent spinup_start (%d) and spinup_stop (%d) values\n",
                   spinup_start_step,
                   spinup_stop_step);
    }
    if (spindown_start_step > -1 || spindown_stop_step > -1) {
        if (spindown_stop_step < spindown_start_step)
            printf("DGRT: inconsistent spindown_start (%d) and spindown_stop (%d) values\n",
                   spindown_start_step,
                   spindown_stop_step);
    }
    RTSetup(Tstar_config,
            planet_star_dist_config,
            radius_star_config,
            diff_ang_config,
            sim.P_Ref,
            sim.Gravit,
            albedo_config,
            esp.kappa_sw,
            esp.kappa_lw,
            latf_lw_config,
            kappa_lw_pole_config,
            n_lw_config,
            n_sw_config,
            esp.f_lw,
            rt1Dmode_config,
            sim.Tmean);

    cudaMemset(surf_flux_d, 0, sizeof(double) * esp.point_num);

    bool returnstatus = true;
    // int  id;
    // if (esp.surface == true) {
    //     if (s != nullptr) {
    //         // load initialisation data from storage s
    //         returnstatus &= (*s).read_table_to_ptr("/Tsurface", Tsurface_h, esp.point_num);
    //     }
    //     else {
    //         for (id = 0; id < esp.point_num; id++) {
    //             Tsurface_h[id] = sim.Tmean;
    //         }
    //         cudaMemset(surf_flux_d, 0, sizeof(double) * esp.point_num);
    //     }
    //     cudaMemcpy(Tsurface_d, Tsurface_h, esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
    // }

    // ask for insolation computation
    esp.insolation.set_require();

    return returnstatus;
}


bool radiative_transfer::phy_loop(ESP &                  esp,
                                  const SimulationSetup &sim,
                                  int                    nstep, // Step number
                                  double                 time_step) {           // Time-step [s]

    bool run      = true;
    Qheat_scaling = 1.0;


    if (spinup_start_step > -1 && spinup_stop_step > -1) {
        if (nstep < spinup_start_step) // before spinup
        {
            run           = false;
            Qheat_scaling = 0.0;
        }
        else if ((nstep >= spinup_start_step) && (nstep <= spinup_stop_step)) // during spinup
        {
            double x = (double)(nstep - spinup_start_step)
                       / (double)(spinup_stop_step - spinup_start_step);
            Qheat_scaling = (1 + sin(M_PI * x - M_PI / 2.0)) / 2.0;
            run           = true;
        }
    }

    if (spindown_start_step > -1 && spindown_stop_step > -1) {
        if ((nstep >= spindown_start_step) && (nstep <= spindown_stop_step)) {
            double x = (double)(nstep - spindown_start_step)
                       / (double)(spindown_stop_step - spindown_start_step);
            Qheat_scaling = 1.0 - (1 + sin(M_PI * x - M_PI / 2.0)) / 2.0;
            run           = true;
        }
        else if (nstep >= spindown_stop_step) {
            run           = false;
            Qheat_scaling = 0.0;
        }
    }

    if (run) {
        //
        //  Number of threads per block.
        const int NTH = 256;

        //  Specify the block sizes.
        dim3 NB((esp.point_num / NTH) + 1, esp.nv, 1);
        dim3 NBRT((esp.point_num / NTH) + 1, 1, 1);

        rtm_dual_band<<<NBRT, NTH>>>(esp.pressure_d,
                                     esp.Rho_d,
                                     esp.temperature_d,
                                     flw_up_d,
                                     flw_dn_d,
                                     fsw_up_d,
                                     fsw_dn_d,
                                     tau_d,
                                     sim.Gravit,
                                     esp.Cp_d,
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
                                     diff_ang,
                                     esp.Tint,
                                     albedo,
                                     tausw,
                                     taulw,
                                     latf_lw,
                                     taulw_pole,
                                     n_sw,
                                     n_lw,
                                     esp.f_lw,
                                     incflx,
                                     sim.P_Ref,
                                     esp.point_num,
                                     esp.nv,
                                     esp.nvi,
                                     sim.A,
                                     esp.insolation.get_r_orb(),
                                     esp.insolation.get_device_cos_zenith_angles(),
                                     insol_d,
                                     esp.surface,
                                     esp.Csurf,
                                     esp.Tsurface_d,
                                     surf_flux_d,
                                     esp.profx_Qheat_d,
                                     qheat_d,
                                     esp.Rd_d,
                                     Qheat_scaling,
                                     sim.gcm_off,
                                     rt1Dmode);

        if (nstep * time_step < (2 * M_PI / esp.insolation.get_mean_motion())) {
            // stationary orbit/obliquity
            // calculate annually average of insolation for the first orbit
            annual_insol<<<NBRT, NTH>>>(insol_ann_d, insol_d, nstep, esp.point_num);
        }
    }
    return true;
}

bool radiative_transfer::configure(config_file &config_reader) {
    // basic star-planet properties
    config_reader.append_config_var("Tstar", Tstar_config, Tstar_config);
    config_reader.append_config_var(
        "planet_star_dist", planet_star_dist_config, planet_star_dist_config);
    config_reader.append_config_var("radius_star", radius_star_config, radius_star_config);
    config_reader.append_config_var("diff_ang", diff_ang_config, diff_ang_config);
    // config_reader.append_config_var("Tint", Tint_config, Tint_config);
    config_reader.append_config_var("albedo", albedo_config, albedo_config);
    // config_reader.append_config_var("tausw", tausw_config, tausw_config);
    // config_reader.append_config_var("taulw", taulw_config, taulw_config);

    // options for latitude dependence in longwave opacity
    config_reader.append_config_var("latf_lw", latf_lw_config, latf_lw_config);
    config_reader.append_config_var("kappa_lw_pole", kappa_lw_pole_config, kappa_lw_pole_config);

    config_reader.append_config_var("n_sw", n_sw_config, n_sw_config);
    config_reader.append_config_var("n_lw", n_lw_config, n_lw_config);
    // config_reader.append_config_var("f_lw", f_lw_config, f_lw_config);


    config_reader.append_config_var("rt1Dmode", rt1Dmode_config, rt1Dmode_config);

    // spin up spin down
    config_reader.append_config_var("dgrt_spinup_start", spinup_start_step, spinup_start_step);
    config_reader.append_config_var("dgrt_spinup_stop", spinup_stop_step, spinup_stop_step);
    config_reader.append_config_var(
        "dgrt_spindown_start", spindown_start_step, spindown_start_step);
    config_reader.append_config_var("dgrt_spindown_stop", spindown_stop_step, spindown_stop_step);

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
        flw_up_h, flw_up_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(flw_up_h, esp.nvi * esp.point_num, "/flw_up", "W m^-2", "upward flux (LW)");

    cudaMemcpy(
        fsw_up_h, fsw_up_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(fsw_up_h, esp.nvi * esp.point_num, "/fsw_up", "W m^-2", "upward flux (SW)");

    cudaMemcpy(
        flw_dn_h, flw_dn_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(flw_dn_h, esp.nvi * esp.point_num, "/flw_dn", "W m^-2", "downward flux (LW)");

    cudaMemcpy(
        fsw_dn_h, fsw_dn_d, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(fsw_dn_h, esp.nvi * esp.point_num, "/fsw_dn", "W m^-2", "downward flux (SW)");

    cudaMemcpy(tau_h, tau_d, esp.nv * esp.point_num * 2 * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(tau_h,
                   esp.nv * esp.point_num * 2,
                   "/tau",
                   " ",
                   "optical depth across each layer (not total optical depth)");

    cudaMemcpy(qheat_h, qheat_d, esp.nv * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table(qheat_h, esp.nv * esp.point_num, "/DGQheat", " ", "Double Gray Qheat");

    // cudaMemcpy(Tsurface_h, Tsurface_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    // s.append_table(Tsurface_h, esp.point_num, "/Tsurface", "K", "surface temperature");

    s.append_value(Qheat_scaling, "/dgrt_qheat_scaling", " ", "Qheat scaling applied to DG");
    return true;
}

bool radiative_transfer::store_init(storage &s) {
    s.append_value(Tstar, "/Tstar", "K", "Temperature of host star");
    // s.append_value(Tint, "/Tint", "K", "Temperature of interior heat flux");
    s.append_value(planet_star_dist / 149597870.7,
                   "/planet_star_dist",
                   "au",
                   "distance b/w host star and planet");
    s.append_value(radius_star / 695508, "/radius_star", "R_sun", "radius of host star");
    s.append_value(diff_ang, "/diff_ang", "-", "diffusivity factor");
    s.append_value(albedo, "/albedo", "-", "bond albedo of planet");
    s.append_value(tausw, "/tausw", "-", "shortwave optical depth of deepest layer");
    s.append_value(taulw, "/taulw", "-", "longwave optical depth of deepest layer");

    s.append_value(latf_lw ? 1.0 : 0.0, "/latf_lw", "-", "use lat dependent opacity");
    s.append_value(
        taulw_pole, "/taulw_pole", "-", "longwave optical depth of deepest layer at poles");
    s.append_value(n_lw, "/n_lw", "-", "power law exponent for unmixed absorber in LW");
    s.append_value(n_sw, "/n_sw", "-", "power law exponent for mixed/unmixed absorber in SW");
    // s.append_value(f_lw, "/f_lw", "-", "fraction of taulw in well-mixed absorber");


    return true;
}

void radiative_transfer::RTSetup(double Tstar_,
                                 double planet_star_dist_,
                                 double radius_star_,
                                 double diff_ang_,
                                 double P_Ref,
                                 double Gravit,
                                 double albedo_,
                                 double kappa_sw,
                                 double kappa_lw,
                                 bool   latf_lw_,
                                 double kappa_lw_pole,
                                 double n_lw_,
                                 double n_sw_,
                                 double f_lw,
                                 bool   rt1Dmode_,
                                 double Tmean) {

    double bc = 5.677036E-8; // Stefan–Boltzmann constant [W m−2 K−4]

    Tstar            = Tstar_;
    planet_star_dist = planet_star_dist_ * 149597870.7; //conv to km
    radius_star      = radius_star_ * 695508;           //conv to km
    diff_ang         = diff_ang_;
    // Tint             = Tint_;
    albedo     = albedo_;
    tausw      = kappa_sw * P_Ref / Gravit;
    taulw      = kappa_lw * P_Ref / (f_lw * Gravit);
    taulw_pole = kappa_lw_pole * P_Ref / (f_lw * Gravit);
    latf_lw    = latf_lw_;
    n_sw       = n_sw_;
    n_lw       = n_lw_;
    // f_lw             = f_lw_;

    double resc_flx = pow(radius_star / planet_star_dist, 2.0);
    incflx          = resc_flx * bc * Tstar * Tstar * Tstar * Tstar;

    rt1Dmode = rt1Dmode_;
}
