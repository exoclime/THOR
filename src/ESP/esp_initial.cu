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
// Build the class ESP (Exoclimes Simulation Platform)
//
//
// Description:
//   Declare and initialize variables in the model
//
// Method: -
//
//
// Known limitations: None.
//
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

#include "directories.h"
#include "esp.h"
#include "log_writer.h"
#include "phy/profx_conservation.h"
#include "phy/ultrahot_thermo.h"
#include "phy/valkyrie_jet_steadystate.h"
#include "storage.h"

#include <map>
#include <stdio.h>

// physical modules
#include "phy_modules.h"

__host__ ESP::ESP(int *                 point_local_,
                  int *                 maps_,
                  double *              lonlat_,
                  double *              Altitude_,
                  double *              Altitudeh_,
                  double *              nvecoa_,
                  double *              nvecti_,
                  double *              nvecte_,
                  double *              areasT_,
                  double *              areasTr_,
                  double *              areas_,
                  double *              div_,
                  double *              grad_,
                  double *              curlz_,
                  double *              func_r_,
                  int                   nl_region_,
                  int                   nr_,
                  int                   nv_,
                  int                   nvi_,
                  int                   glevel_,
                  bool                  spring_dynamics_,
                  double                spring_beta_,
                  int                   nlat_bins_,
                  int *                 zonal_mean_tab,
                  double                Ruv_sponge_,
                  double                Rw_sponge_,
                  double                RT_sponge_,
                  double                ns_ray_sponge_,
                  bool                  damp_uv_to_mean_,
                  bool                  damp_w_to_mean_,
                  raysp_calc_mode_types raysp_calc_mode_,
                  double                Dv_sponge_,
                  double                ns_diff_sponge_,
                  int                   order_diff_sponge_,
                  double                t_shrink_,
                  int                   point_num_,
                  bool                  conservation,
                  benchmark_types       core_benchmark_,
                  log_writer &          logwriter_,
                  int                   max_count_,
                  bool                  output_mean,
                  init_PT_profile_types init_PT_profile_,
                  double                Tint_,
                  double                kappa_lw_,
                  double                kappa_sw_,
                  double                f_lw_,
                  uh_thermo_types       ultrahot_thermo_,
                  uh_heating_types      ultrahot_heating_,
                  thermo_equation_types thermo_equation_) :
    nl_region(nl_region_),
    nr(nr_),
    point_num(point_num_),
    nv(nv_),
    nvi(nvi_),
    nlat_bins(nlat_bins_),
    order_diff_sponge(order_diff_sponge_),
    damp_uv_to_mean(damp_uv_to_mean_),
    damp_w_to_mean(damp_w_to_mean_),
    glevel(glevel_),
    spring_dynamics(spring_dynamics_),
    spring_beta(spring_beta_),
    logwriter(logwriter_),
    core_benchmark(core_benchmark_),
    init_PT_profile(init_PT_profile_),
    raysp_calc_mode(raysp_calc_mode_),
    ultrahot_thermo(ultrahot_thermo_),
    ultrahot_heating(ultrahot_heating_),
    thermo_equation(thermo_equation_) {

    point_local_h = point_local_;
    maps_h        = maps_;

    lonlat_h = lonlat_;

    Altitude_h  = Altitude_;
    Altitudeh_h = Altitudeh_;

    nvecoa_h  = nvecoa_;
    nvecti_h  = nvecti_;
    nvecte_h  = nvecte_;
    areasTr_h = areasTr_;
    areasT_h  = areasT_;
    areas_h   = areas_;

    div_h   = div_;
    grad_h  = grad_;
    curlz_h = curlz_;

    func_r_h = func_r_;

    zonal_mean_tab_h = zonal_mean_tab;

    Ruv_sponge     = Ruv_sponge_;
    Rw_sponge      = Rw_sponge_;
    RT_sponge      = RT_sponge_;
    ns_ray_sponge  = ns_ray_sponge_;
    Dv_sponge      = Dv_sponge_;
    ns_diff_sponge = ns_diff_sponge_;

    t_shrink  = t_shrink_;
    max_count = max_count_;

    Tint     = Tint_;
    kappa_lw = kappa_lw_;
    kappa_sw = kappa_sw_;
    f_lw     = f_lw_;

    // Set the physics module execute state for the rest of the lifetime of ESP object
    // only execute physics modules when no benchmarks are enabled
    if (core_benchmark == NO_BENCHMARK) {
        phy_modules_execute = true;
    }
    else
        phy_modules_execute = false;

    //
    //  Allocate Data
    alloc_data(conservation, output_mean);
}

__host__ void ESP::alloc_data(bool conservation, bool output_mean) {

    //
    //  Description:
    //
    //  Allocate data on host and device.
    //
    //  Allocate data in host
    //  Diagnostics an doutput
    Rho_h         = (double *)malloc(nv * point_num * sizeof(double));
    pressure_h    = (double *)malloc(nv * point_num * sizeof(double));
    temperature_h = (double *)malloc(nv * point_num * sizeof(double));
    Mh_h          = (double *)malloc(nv * point_num * 3 * sizeof(double));
    W_h           = (double *)malloc(nv * point_num * sizeof(double));
    Wh_h          = (double *)malloc(nvi * point_num * sizeof(double));

    if (output_mean == true) {
        Rho_mean_h      = (double *)malloc(nv * point_num * sizeof(double));
        pressure_mean_h = (double *)malloc(nv * point_num * sizeof(double));
        Mh_mean_h       = (double *)malloc(nv * point_num * 3 * sizeof(double));
        Wh_mean_h       = (double *)malloc(nvi * point_num * sizeof(double));
    }

    if (conservation == true) {
        Etotal_h  = (double *)malloc(nv * point_num * sizeof(double));
        Mass_h    = (double *)malloc(nv * point_num * sizeof(double));
        AngMomx_h = (double *)malloc(nv * point_num * sizeof(double));
        AngMomy_h = (double *)malloc(nv * point_num * sizeof(double));
        AngMomz_h = (double *)malloc(nv * point_num * sizeof(double));
        Entropy_h = (double *)malloc(nv * point_num * sizeof(double));
    }

    // ultra-hot jupiter stuff
    Rd_h = (double *)malloc(nv * point_num * sizeof(double));
    Cp_h = (double *)malloc(nv * point_num * sizeof(double));

    flux_vec        = (double *)malloc(nv * point_num * sizeof(double));
    boundary_flux_h = (double *)malloc(6 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&boundary_flux_d, 6 * point_num * nv * sizeof(double));

    //  Allocate data in device
    //  Grid
    cudaMalloc((void **)&point_local_d, 6 * point_num * sizeof(int));
    cudaMalloc((void **)&maps_d, (nl_region + 2) * (nl_region + 2) * nr * sizeof(int));

    //  Operators
    cudaMalloc((void **)&nvecoa_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecti_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecte_d, 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&areasT_d, point_num * sizeof(double));
    cudaMalloc((void **)&areasTr_d, 6 * point_num * sizeof(double));
    cudaMalloc((void **)&areas_d, 3 * 6 * point_num * sizeof(double));
    cudaMalloc((void **)&func_r_d, 3 * point_num * sizeof(double));
    cudaMalloc((void **)&div_d, 7 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&grad_d, 7 * 3 * point_num * sizeof(double));

    //  Altitude (grid)
    cudaMalloc((void **)&Altitude_d, nv * sizeof(double));
    cudaMalloc((void **)&Altitudeh_d, nvi * sizeof(double));

    //  Longitude-latitude
    cudaMalloc((void **)&lonlat_d, 2 * point_num * sizeof(double));

    //  Diagnostics
    cudaMalloc((void **)&Mh_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&W_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Wh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&Rho_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&pressure_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&pressureh_d, (nv + 1) * point_num * sizeof(double));

    if (output_mean == true) {
        // Average quantities over output interval
        cudaMalloc((void **)&Mh_mean_d, nv * point_num * 3 * sizeof(double));
        cudaMalloc((void **)&Wh_mean_d, nvi * point_num * sizeof(double));
        cudaMalloc((void **)&Rho_mean_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&pressure_mean_d, nv * point_num * sizeof(double));
    }

    // ultra hot
    cudaMalloc((void **)&Rd_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Cp_d, nv * point_num * sizeof(double));

    //  Temperature
    cudaMalloc((void **)&temperature_d, nv * point_num * sizeof(double));

    //  Potential temperature
    cudaMalloc((void **)&pt_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&pth_d, nvi * point_num * sizeof(double));

    //  Energy (for thermo_equation = energy)
    cudaMalloc((void **)&epotential_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&epotentialh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&ekinetic_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&ekinetich_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&Etotal_tau_d, nv * point_num * sizeof(double));

    //  Entalphy
    cudaMalloc((void **)&h_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&hh_d, nvi * point_num * sizeof(double));

    //  Advection
    cudaMalloc((void **)&Adv_d, nv * point_num * 3 * sizeof(double));

    //  3D vector
    cudaMalloc((void **)&v_d, nv * point_num * 3 * sizeof(double));

    //  Effective gravity
    cudaMalloc((void **)&gtil_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&gtilh_d, nvi * point_num * sizeof(double));

    //  Slow modes
    cudaMalloc((void **)&SlowMh_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&SlowWh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&SlowRho_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Slowpressure_d, nv * point_num * sizeof(double));


    //  Deviations
    cudaMalloc((void **)&pressures_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Rhos_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Mhs_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Ws_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Whs_d, nvi * point_num * sizeof(double));


    //  RK-Method
    cudaMalloc((void **)&pressurek_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Rhok_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Mhk_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Wk_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Whk_d, nvi * point_num * sizeof(double));

    //  Vertical integration
    cudaMalloc((void **)&Sp_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&Sd_d, nv * point_num * sizeof(double));

    //  Diffusion
    cudaMalloc((void **)&Kdhz_d, nv * sizeof(double));
    cudaMalloc((void **)&Kdh4_d, nv * sizeof(double));
    cudaMalloc((void **)&Kdvz_d, nv * sizeof(double));
    cudaMalloc((void **)&Kdv6_d, nv * sizeof(double));

    cudaMalloc((void **)&DivM_d, nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&diffpr_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffmh_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffw_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffrh_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diff_d, 6 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&divg_Mh_d, 3 * nv * point_num * sizeof(double));

    cudaMalloc((void **)&Kdh2_d, nv * sizeof(double));

    cudaMalloc((void **)&diffprv_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffmv_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffwv_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffrv_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffv_d1, 6 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&diffv_d2, 6 * nv * point_num * sizeof(double));

    cudaMalloc((void **)&profx_Qheat_d, nv * point_num * sizeof(double));
    cudaMalloc((void **)&profx_dMh_d, 3 * nv * point_num * sizeof(double));
    cudaMalloc((void **)&profx_dWh_d, nvi * point_num * sizeof(double));
    cudaMalloc((void **)&profx_dW_d, nv * point_num * sizeof(double));

    //  Extras-nan
    cudaMalloc((void **)&check_d, sizeof(bool));

    cudaMalloc((void **)&vbar_d, 3 * nv * nlat_bins * sizeof(double));
    cudaMalloc((void **)&zonal_mean_tab_d, 3 * point_num * sizeof(int));
    vbar_h = (double *)malloc(3 * nv * nlat_bins * sizeof(double));
    cudaMalloc((void **)&utmp, nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&vtmp, nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&wtmp, nv * nlat_bins * max_count * sizeof(double));
    utmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));
    vtmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));
    wtmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));
    cudaMalloc((void **)&Tbar_d, nv * nlat_bins * sizeof(double));
    Tbar_h = (double *)malloc(nv * nlat_bins * sizeof(double));
    cudaMalloc((void **)&Ttmp, nv * nlat_bins * max_count * sizeof(double));
    Ttmp_h = (double *)malloc(nv * nlat_bins * max_count * sizeof(double));

    if (conservation == true) {
        //  Conservation quantities
        cudaMalloc((void **)&Etotal_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&Entropy_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&Mass_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomx_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomy_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&AngMomz_d, nv * point_num * sizeof(double));
        cudaMalloc((void **)&GlobalE_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalEnt_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalMass_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMx_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMy_d, 1 * sizeof(double));
        cudaMalloc((void **)&GlobalAMz_d, 1 * sizeof(double));
    }

    // PHY modules
    log::printf("  Dynamical core memory initialised.\n");

    if (phy_modules_execute) {

        // physics module need to initialise their own memory
        bool init_modules = phy_modules_init_mem(*this, phy_modules_core_arrays);
        // Physics module register arrays that need to be updated in dynamical core Runge-Kutta step
        phy_modules_core_arrays.allocate_device_array();
        if (init_modules)
            log::printf("  Module memory initialised.\n");
        else {
            log::printf("  Error initialising module memory.\n");
            exit(-1);
        }
    }
}

__host__ bool ESP::initial_values(const std::string &initial_conditions_filename,
                                  const bool &       continue_sim,
                                  double             timestep_dyn,
                                  SimulationSetup &  sim,
                                  int &              nstep,
                                  double &           simulation_start_time,
                                  int &              output_file_idx) {

    output_file_idx = 0;
    nstep           = 0;
    //  Set initial conditions.
    //
    //
    //  Initial atmospheric conditions
    bool   read_gibbs = read_in_gibbs_H(GibbsN); //ultrahot jup
    double chi_H = 0, ptmp, eps = 1e-8, f, df, dz, mu;
    int    it, it_max = 100;

    double Rd_L, P_L, T_L;
    if (sim.rest) {
        for (int i = 0; i < point_num; i++) {
            //
            //          Initial conditions for an isothermal Atmosphere
            //
            if (init_PT_profile == ISOTHERMAL && ultrahot_thermo == NO_UH_THERMO) {
                //isothermal initial profile, no variation in Rd or Cp due to H-H2 reaction
                double Ha = sim.Rd * sim.Tmean / sim.Gravit;
                for (int lev = 0; lev < nv; lev++) {
                    pressure_h[i * nv + lev]    = sim.P_Ref * exp(-Altitude_h[lev] / Ha);
                    temperature_h[i * nv + lev] = sim.Tmean;
                    Rd_h[i * nv + lev]          = sim.Rd;
                    Cp_h[i * nv + lev]          = sim.Cp;
                }
            }
            else {
                //
                //          Initial conditions for a non-isothermal Atmosphere
                //
                mu = 0.5;

                for (int lev = 0; lev < nv; lev++) {
                    if (lev == 0) {
                        if (init_PT_profile == ISOTHERMAL) {
                            temperature_h[i * nv + lev] = sim.Tmean;
                        }
                        else {
                            temperature_h[i * nv + lev] = guillot_T(sim.P_Ref,
                                                                    mu,
                                                                    sim.Tmean,
                                                                    sim.P_Ref,
                                                                    sim.Gravit,
                                                                    Tint,
                                                                    f_lw,
                                                                    kappa_sw,
                                                                    kappa_lw);
                        }
                        if (ultrahot_thermo != NO_UH_THERMO) {
                            chi_H = chi_H_equilibrium(
                                GibbsT, GibbsdG, GibbsN, temperature_h[i * nv + lev], sim.P_Ref);
                            Rd_L = Rd_from_chi_H(chi_H);
                        }
                        else {
                            Rd_L = sim.Rd;
                        }
                        P_L = sim.P_Ref;
                        T_L = temperature_h[i * nv + lev];
                        dz  = Altitude_h[0];
                    }
                    else {
                        temperature_h[i * nv + lev] = temperature_h[i * nv + lev - 1];
                        if (ultrahot_thermo != NO_UH_THERMO) {
                            chi_H = chi_H_equilibrium(
                                GibbsT, GibbsdG, GibbsN, sim.Tmean, pressure_h[i * nv + lev - 1]);
                            Rd_L = Rd_h[i * nv + lev - 1];
                        }
                        else {
                            Rd_L = Rd_h[i * nv + lev - 1];
                        }
                        P_L = pressure_h[i * nv + lev - 1];
                        T_L = temperature_h[i * nv + lev - 1];
                        dz  = Altitude_h[lev] - Altitude_h[lev - 1];
                    }
                    pressure_h[i * nv + lev] = P_L;
                    Rd_h[i * nv + lev]       = Rd_L;
                    ptmp                     = pressure_h[i * nv + lev] + 2 * eps;

                    it = 0;
                    while (it < it_max && ptmp - pressure_h[i * nv + lev] > eps) {
                        ptmp = pressure_h[i * nv + lev];
                        f    = log(pressure_h[i * nv + lev] / P_L) / dz
                            + sim.Gravit
                                  / (0.5
                                     * (Rd_h[i * nv + lev] * temperature_h[i * nv + lev]
                                        + Rd_L * T_L));
                        df                       = 1.0 / (pressure_h[i * nv + lev] * dz);
                        pressure_h[i * nv + lev] = pressure_h[i * nv + lev] - f / df;
                        if (init_PT_profile == ISOTHERMAL) {
                            temperature_h[i * nv + lev] = sim.Tmean;
                        }
                        else {
                            temperature_h[i * nv + lev] = guillot_T(pressure_h[i * nv + lev],
                                                                    mu,
                                                                    sim.Tmean,
                                                                    sim.P_Ref,
                                                                    sim.Gravit,
                                                                    Tint,
                                                                    f_lw,
                                                                    kappa_sw,
                                                                    kappa_lw);
                        }
                        if (ultrahot_thermo != NO_UH_THERMO) {
                            chi_H              = chi_H_equilibrium(GibbsT,
                                                      GibbsdG,
                                                      GibbsN,
                                                      temperature_h[i * nv + lev],
                                                      pressure_h[i * nv + lev]);
                            Rd_h[i * nv + lev] = Rd_from_chi_H(chi_H);
                        }
                        else {
                            Rd_h[i * nv + lev] = sim.Rd;
                        }
                        it++;
                    }
                    if (ultrahot_thermo != NO_UH_THERMO) {
                        Cp_h[i * nv + lev] = Cp_from_chi_H(chi_H, temperature_h[i * nv + lev]);
                    }
                    else {
                        Cp_h[i * nv + lev] = sim.Cp;
                    }
                }
            }

            for (int lev = 0; lev < nv; lev++) {
                //              Density [kg/m3]
                Rho_h[i * nv + lev] =
                    pressure_h[i * nv + lev] / (temperature_h[i * nv + lev] * Rd_h[i * nv + lev]);

                //              Momentum [kg/m3 m/s]
                Mh_h[i * 3 * nv + 3 * lev + 0] = 0.0;
                Mh_h[i * 3 * nv + 3 * lev + 1] = 0.0;
                Mh_h[i * 3 * nv + 3 * lev + 2] = 0.0;

                //              Vertical momentum [kg/m3 m/s]
                W_h[i * nv + lev]        = 0.0; // Center of the layer.
                Wh_h[i * (nv + 1) + lev] = 0.0; // Layers interface.
            }
            Wh_h[i * (nv + 1) + nv] = 0.0;
        }
        if (core_benchmark == JET_STEADY) {
            //  Number of threads per block.
            const int NTH = 256;

            //  Specify the block sizes.
            dim3 NB((point_num / NTH) + 1, nv, 1);

            cudaMemcpy(Altitude_d, Altitude_h, nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(
                pressure_d, pressure_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(Mh_d, Mh_h, 3 * point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(Rho_d, Rho_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(temperature_d,
                       temperature_h,
                       point_num * nv * sizeof(double),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(lonlat_d, lonlat_h, 2 * point_num * sizeof(double), cudaMemcpyHostToDevice);
            setup_jet<<<NB, NTH>>>(Mh_d,
                                   // setup_jet <<< 1, 1 >>>  (Mh_d,
                                   pressure_d,
                                   Rho_d,
                                   temperature_d,
                                   sim.Cp,
                                   sim.Rd,
                                   sim.Omega,
                                   sim.A,
                                   Altitude_d,
                                   lonlat_d,
                                   point_num);

            cudaMemcpy(Mh_h, Mh_d, 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(temperature_h,
                       temperature_d,
                       point_num * nv * sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(
                pressure_h, pressure_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(Rho_h, Rho_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        }

        simulation_start_time = 0.0;
    } //end if rest
    else {
        bool load_OK = true;
        // build planet filename
        string planet_filename;

        path   p(initial_conditions_filename);
        int    file_number = 0;
        string basename    = "";

        string parent_path = p.parent();

        // Reload correct file if we are continuing from a specific file
        if (continue_sim) {
            if (!match_output_file_numbering_scheme(
                    initial_conditions_filename, basename, file_number)) {
                log::printf("Loading initial conditions: "
                            "Could not recognise file numbering scheme "
                            "for input %s: (found base: %s, num: %d) \n",
                            initial_conditions_filename.c_str(),
                            basename.c_str(),
                            file_number);
                return false;
            }

            output_file_idx = file_number;

            planet_filename = p.parent() + "/esp_output_planet_" + basename + ".h5";
        }
        else {
            planet_filename = p.parent() + "/" + p.stem() + "_planet.h5";
        }

        // check existence of files
        if (!path_exists(initial_conditions_filename)) {
            log::printf("initial condition file %s not found.\n",
                        initial_conditions_filename.c_str());
            return false;
        }

        if (!path_exists(planet_filename)) {
            log::printf("planet_file %s not found.\n", planet_filename.c_str());
            return false;
        }


        log::printf("Loading planet from: %s\n", planet_filename.c_str());
        log::printf("Loading initial conditions from: %s\n", initial_conditions_filename.c_str());

        // Check planet data
        {
            // values from initial conditions to check against variables from config
            map<string, double> mapValuesDouble;
            map<string, int>    mapValuesInt;

            mapValuesDouble["/A"]            = sim.A;
            mapValuesDouble["/Top_altitude"] = sim.Top_altitude;
            mapValuesInt["/glevel"]          = glevel;
            mapValuesInt["/vlevel"]          = nv;

            storage s(planet_filename, true);

            bool values_match = true;
            // double

            for (const std::pair<std::string, double> &element : mapValuesDouble) {
                double value = 0.0;
                load_OK      = s.read_value(element.first, value);

                if (!load_OK) {
                    printf("Error reading key %s from reload config.\n", element.first.c_str());
                    values_match = false;
                }


                if (value != element.second) {
                    log::printf("mismatch for %s value between config value: %f and initial "
                                "condition value %f.\n",
                                element.first.c_str(),
                                element.second,
                                value);
                    values_match = false;
                }
            }

            // int var
            for (const std::pair<std::string, int> &element : mapValuesInt) {
                int value = 0;
                load_OK   = s.read_value(element.first, value);

                if (!load_OK) {
                    printf("Error reading key %s from reload config.\n", element.first.c_str());
                    values_match = false;
                }

                if (value != element.second) {
                    log::printf("mismatch for %s value between config value: %d and initial "
                                "condition value %d.\n",
                                element.first.c_str(),
                                element.second,
                                value);
                    values_match = false;
                }
            }


            if (load_OK == false || values_match == false) {
                log::printf("Could not reload full configuration.\n");

                return false;
            }
        }


        //      Restart from an existing simulation.
        {
            // Load atmospheric data
            storage s(initial_conditions_filename, true);
            // Step number
            load_OK &= s.read_value("/nstep", nstep);

            log::printf("Reloaded %s: %d.\n", "/nstep", load_OK ? 1 : 0);

            //      Density
            load_OK &= s.read_table_to_ptr("/Rho", Rho_h, point_num * nv);
            log::printf("Reloaded %s: %d.\n", "/Rho", load_OK ? 1 : 0);
            //      Pressure
            load_OK &= s.read_table_to_ptr("/Pressure", pressure_h, point_num * nv);
            log::printf("Reloaded %s: %d.\n", "/Pressure", load_OK ? 1 : 0);
            //      Horizontal momentum
            load_OK &= s.read_table_to_ptr("/Mh", Mh_h, point_num * nv * 3);
            log::printf("Reloaded %s: %d.\n", "/Mh", load_OK ? 1 : 0);
            //      Vertical momentum
            load_OK &= s.read_table_to_ptr("/Wh", Wh_h, point_num * nvi);
            log::printf("Reloaded %s: %d.\n", "/Wh", load_OK ? 1 : 0);

            load_OK &= s.read_table_to_ptr("/Rd", Rd_h, point_num * nv);
            log::printf("Reloaded %s: %d.\n", "/Rd", load_OK ? 1 : 0);

            load_OK &= s.read_table_to_ptr("/Cp", Cp_h, point_num * nv);
            log::printf("Reloaded %s: %d.\n", "/Cp", load_OK ? 1 : 0);
            //      Simulation start time
            load_OK &= s.read_value("/simulation_time", simulation_start_time);
            log::printf("Reloaded %s: %d.\n", "/simulation_time", load_OK ? 1 : 0);
        }


        if (!load_OK) {
            log::printf("Error reloading simulation state\n");

            return false;
        }


        for (int i = 0; i < point_num; i++) {
            for (int lev = 0; lev < nv; lev++) {
                //hack
                // Rd_h[i * nv + lev] = sim.Rd;
                // Cp_h[i * nv + lev] = sim.Cp;
                ////
                temperature_h[i * nv + lev] =
                    pressure_h[i * nv + lev] / (Rd_h[i * nv + lev] * Rho_h[i * nv + lev]);
            }
        }


        for (int i = 0; i < point_num; i++) {
            for (int lev = 0; lev < nv; lev++) {
                double xi   = Altitude_h[lev];
                double xim1 = Altitudeh_h[lev];
                double xip1 = Altitudeh_h[lev + 1];

                double a = (xi - xip1) / (xim1 - xip1);
                double b = (xi - xim1) / (xip1 - xim1);

                W_h[i * nv + lev] = Wh_h[i * (nv + 1) + lev] * a + Wh_h[i * (nv + 1) + lev + 1] * b;
            }
        }
    } //end if rest == false
#ifdef BENCHMARKING
    // recompute temperature from pressure and density, to have correct rounding for binary comparison
    for (int i = 0; i < point_num; i++)
        for (int lev = 0; lev < nv; lev++)
            temperature_h[i * nv + lev] =
                pressure_h[i * nv + lev] / (Rd_h[i * nv + lev] * Rho_h[i * nv + lev]);
#endif // BENCHMARKING

    //  Diffusion
    //  Horizontal
    double *Kdhz_h, *Kdh4_h;
    Kdhz_h = new double[nv]; // horizontal divergence damping strength
    Kdh4_h = new double[nv]; // horizontal diffusion strength
                             // if (sim.DiffSponge) {
    double  n, ksponge;
    double *Kdh2_h;
    Kdh2_h = new double[nv];
    for (int lev = 0; lev < nv; lev++) {
        double dbar = sqrt(2 * M_PI / 5) * sim.A / (pow(2, glevel));
        Kdh4_h[lev] =
            (sim.Diffc) * pow(dbar, 4.) / timestep_dyn; // * Altitude_h[lev]/sim.Top_altitude;
        Kdhz_h[lev] =
            (sim.DivDampc) * pow(dbar, 4.) / timestep_dyn; // * Altitude_h[lev]/sim.Top_altitude;
        if (sim.DiffSponge) {
            n = Altitude_h[lev] / sim.Top_altitude;
            if (n > ns_diff_sponge) {
                ksponge = Dv_sponge
                          * pow(sin(0.5 * M_PI * (n - ns_diff_sponge) / (1.0 - ns_diff_sponge)), 2);
            }
            else {
                ksponge = 0;
            }
            if (order_diff_sponge == 2) {
                Kdh2_h[lev] = ksponge * pow(dbar, 2.) / timestep_dyn;
            }
            else if (order_diff_sponge == 4) {
                Kdh4_h[lev] += ksponge * pow(dbar, 4.) / timestep_dyn;
            }
        }
    }

    //  Diffusion
    //  Vertical
    double *Kdvz_h, *Kdv6_h;
    Kdvz_h = new double[nv]; // vertical divergence damping strength
    Kdv6_h = new double[nv]; // vertical diffusion strength
    for (int lev = 0; lev < nv; lev++) {
        //      Diffusion constant.
        // double dz   = sim.Top_altitude / nv;
        Kdv6_h[lev] = 0.0; //not used (yet? perhaps in future)
        Kdvz_h[lev] = 0.0; //not used (yet? perhaps in future)
    }


    //  Copy memory to the device
    cudaMemcpy(point_local_d, point_local_h, 6 * point_num * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(maps_d,
               maps_h,
               (nl_region + 2) * (nl_region + 2) * nr * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(Altitude_d, Altitude_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Altitudeh_d, Altitudeh_h, nvi * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecoa_d, nvecoa_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecti_d, nvecti_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecte_d, nvecte_h, 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areasTr_d, areasTr_h, 6 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areasT_d, areasT_h, point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areas_d, areas_h, 3 * 6 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(lonlat_d, lonlat_h, 2 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(func_r_d, func_r_h, 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(
        temperature_d, temperature_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Mh_d, Mh_h, point_num * nv * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(W_d, W_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Wh_d, Wh_h, point_num * nvi * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Rho_d, Rho_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pressure_d, pressure_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(div_d, div_h, 7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(grad_d, grad_h, 7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdhz_d, Kdhz_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdh4_d, Kdh4_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdvz_d, Kdvz_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdv6_d, Kdv6_h, nv * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(Kdh2_d, Kdh2_h, nv * sizeof(double), cudaMemcpyHostToDevice);


    if (sim.output_mean == true) {
        cudaMemcpy(Mh_mean_d, Mh_h, point_num * nv * 3 * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(
            pressure_mean_d, pressure_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(Wh_mean_d, Wh_h, point_num * nvi * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(Rho_mean_d, Rho_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    }

    if (sim.RayleighSponge == true)
        cudaMemcpy(zonal_mean_tab_d,
                   zonal_mean_tab_h,
                   3 * point_num * sizeof(int),
                   cudaMemcpyHostToDevice);

    cudaMemcpy(Rd_d, Rd_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Cp_d, Cp_h, point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(GibbsT_d, GibbsT, GibbsN * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(GibbsdG_d, GibbsdG, GibbsN * sizeof(double), cudaMemcpyHostToDevice);

    //  Initialize arrays
    cudaMemset(Adv_d, 0, sizeof(double) * 3 * point_num * nv);
    cudaMemset(v_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(pt_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pth_d, 0, sizeof(double) * nvi * point_num);
    // cudaMemset(pt_tau_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(epotential_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(epotentialh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(ekinetic_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(ekinetich_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Etotal_tau_d, 0, sizeof(double) * nv * point_num);

    cudaMemset(SlowMh_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(SlowWh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(SlowRho_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Slowpressure_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(h_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(hh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Rhos_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pressures_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Mhs_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Ws_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Whs_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(gtil_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(gtilh_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Rhok_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(pressurek_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Mhk_d, 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Wk_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(Whk_d, 0, sizeof(double) * nvi * point_num);
    cudaMemset(Sp_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(Sd_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(DivM_d, 0, sizeof(double) * point_num * 3 * nv);
    cudaMemset(diffpr_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffmh_d, 0, sizeof(double) * 3 * nv * point_num);
    cudaMemset(diffw_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffrh_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diff_d, 0, sizeof(double) * 6 * nv * point_num);
    cudaMemset(divg_Mh_d, 0, sizeof(double) * 3 * nv * point_num);

    cudaMemset(diffprv_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffmv_d, 0, sizeof(double) * 3 * nv * point_num);
    cudaMemset(diffwv_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffrv_d, 0, sizeof(double) * nv * point_num);
    cudaMemset(diffv_d1, 0, sizeof(double) * 6 * nv * point_num);
    cudaMemset(diffv_d2, 0, sizeof(double) * 6 * nv * point_num);

    cudaMemset(boundary_flux_d, 0, sizeof(double) * 6 * nv * point_num);

    delete[] Kdh4_h;
    delete[] Kdhz_h;
    delete[] Kdv6_h;
    delete[] Kdvz_h;
    delete[] Kdh2_h;

    // modules need to set their initial conditions
    if (phy_modules_execute) {
        if (sim.rest) // no initial condition file
            phy_modules_init_data(*this, sim, nullptr);
        else {
            // load initial condition file and pass it to modules
            storage s(initial_conditions_filename, true);

            phy_modules_init_data(*this, sim, &s);
        }
    }


    return true;
}

__host__ ESP::~ESP() {

    //
    //  Description: Frees the memory space.
    //
    //  Host
    // Simulation state data
    free(Rho_h);
    free(pressure_h);
    free(temperature_h);
    free(Mh_h);
    free(W_h);
    free(Wh_h);

    //  Device
    // Grid
    cudaFree(point_local_d);
    cudaFree(maps_d);

    //  Altitude (grid)
    cudaFree(Altitude_d);
    cudaFree(Altitudeh_d);

    //  Operators
    cudaFree(nvecoa_d);
    cudaFree(nvecti_d);
    cudaFree(nvecte_d);
    cudaFree(areasT_d);
    cudaFree(areasTr_d);
    cudaFree(areas_d);

    //  Longitude-latitude
    cudaFree(lonlat_d);
    cudaFree(div_d);
    cudaFree(grad_d);
    cudaFree(func_r_d);

    //  Temperature
    cudaFree(temperature_d);
    //  Diagnostics
    cudaFree(Mh_d);

    cudaFree(W_d);
    cudaFree(Wh_d);
    cudaFree(Rho_d);
    cudaFree(pressure_d);
    cudaFree(pressureh_d);

    //  Entalphy
    cudaFree(h_d);
    cudaFree(hh_d);

    //  Advection
    cudaFree(Adv_d);
    //  Effective gravity
    cudaFree(gtil_d);
    cudaFree(gtilh_d);
    //  3D vector
    cudaFree(v_d);
    //  Potential temperature
    cudaFree(pt_d);
    cudaFree(pth_d);
    //  Slow modes
    cudaFree(SlowMh_d);
    cudaFree(SlowWh_d);
    cudaFree(SlowRho_d);
    cudaFree(Slowpressure_d);
    //  RK-Method
    cudaFree(Rhok_d);
    cudaFree(pressurek_d);
    cudaFree(Mhk_d);
    cudaFree(Whk_d);
    cudaFree(Wk_d);
    //  Deviations
    cudaFree(Rhos_d);
    cudaFree(pressures_d);
    cudaFree(Mhs_d);
    cudaFree(Whs_d);
    cudaFree(Ws_d);

    //  Vertical integration
    cudaFree(Sd_d);
    cudaFree(Sp_d);

    //  Diffusion
    cudaFree(Kdhz_d);
    cudaFree(Kdh4_d);
    cudaFree(Kdvz_d);
    cudaFree(Kdv6_d);
    cudaFree(DivM_d);
    cudaFree(diffpr_d);
    cudaFree(diffmh_d);
    cudaFree(diffw_d);
    cudaFree(diffrh_d);
    cudaFree(diff_d);
    cudaFree(divg_Mh_d);

    cudaFree(Kdh2_d);

    cudaFree(diffprv_d);
    cudaFree(diffmv_d);
    cudaFree(diffwv_d);
    cudaFree(diffrv_d);
    cudaFree(diffv_d1);
    cudaFree(diffv_d2);

    //  Conservation quantities
    cudaFree(Etotal_d);
    cudaFree(Entropy_d);
    cudaFree(Mass_d);
    cudaFree(AngMomx_d);
    cudaFree(AngMomy_d);
    cudaFree(AngMomz_d);
    cudaFree(GlobalE_d);
    cudaFree(GlobalEnt_d);
    cudaFree(GlobalMass_d);
    cudaFree(GlobalAMx_d);
    cudaFree(GlobalAMy_d);
    cudaFree(GlobalAMz_d);
    free(Etotal_h);
    free(Entropy_h);
    free(Mass_h);
    free(AngMomx_h);
    free(AngMomy_h);
    free(AngMomz_h);

    //  Extras-nan
    cudaFree(check_d);

    // Sponge Layer
    cudaFree(vbar_d);
    cudaFree(zonal_mean_tab_d);
    cudaFree(Tbar_d);

    free(vbar_h);
    free(utmp_h);
    free(vtmp_h);
    free(wtmp_h);

    free(Tbar_h);
    free(Ttmp_h);

    cudaFree(utmp);
    cudaFree(vtmp);
    cudaFree(wtmp);
    cudaFree(Ttmp);

    cudaFree(profx_Qheat_d);
    cudaFree(profx_dMh_d);
    cudaFree(profx_dWh_d);
    cudaFree(profx_dW_d);

    cudaFree(epotential_d);
    cudaFree(epotentialh_d);
    cudaFree(ekinetic_d);
    cudaFree(ekinetich_d);
    cudaFree(Etotal_tau_d);

    // ultra hot
    free(Rd_h);
    free(Cp_h);
    free(GibbsT);
    free(GibbsdG);

    cudaFree(Rd_d);
    cudaFree(Cp_d);

    cudaFree(GibbsT_d);
    cudaFree(GibbsdG_d);

    cudaFree(boundary_flux_d);
    free(boundary_flux_h);

    if (phy_modules_execute)
        phy_modules_free_mem();


    log::printf("\n\n Free memory!\n\n");
}
