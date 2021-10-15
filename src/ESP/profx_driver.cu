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
// Description: Physics modules.
//
//
// Method: This version just includes the held-suarez test.
//
// Known limitations: None
//
// Known issues: None
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
#include "phy/dry_conv_adj.h"
#include "phy/profx_acoustic_test.h"
#include "phy/profx_auxiliary.h"
#include "phy/profx_deepHJ.h"
#include "phy/profx_globdiag.h"
#include "phy/profx_gwave_test.h"
#include "phy/profx_held_suarez.h"
#include "phy/profx_shallowHJ.h"
#include "phy/profx_sponge.h"
#include "phy/profx_tidalearth.h"
#include "phy/ultrahot_thermo.h"

#include "binary_test.h"
#include "debug_helpers.h"
#include "log_writer.h"
#include "phy_modules.h"

#include "reduction_add.h"

#include "diagnostics.h"
#include "insolation.h"

__host__ void ESP::ProfX(const SimulationSetup& sim,
                         int                    n_out,
                         kernel_diagnostics&    diag) { // output step (triggers globdiag calc)
    USE_BENCHMARK()
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3      NB((point_num / NTH) + 1, nv, 1);
    dim3      NBRT((point_num / NTH) + 1, 1, 1);

    cudaMemset(profx_Qheat_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(dTsurf_dt_d, 0, sizeof(double) * point_num);

    if (sim.out_interm_momentum) {
        cudaMemcpy(
            Mh_start_dt_d, Mh_d, point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    }

    Recompute_W<<<NB, NTH>>>(
        W_d, Wh_d, Altitude_d, Altitudeh_d, point_num); //trying to stamp out bincomp issue

    cudaDeviceSynchronize();
    if (sim.RayleighSponge == true && raysp_calc_mode == IMP) {
        dim3 NBT((point_num / NTH) + 1, nv, 1);

        if (damp_uv_to_mean) {
            cudaMemset(utmp, 0, sizeof(double) * nlat_bins * nv * max_count);
            cudaMemset(vtmp, 0, sizeof(double) * nlat_bins * nv * max_count);
            zonal_uv<<<NB, NTH>>>(
                Mh_d, Rho_d, zonal_mean_tab_d, lonlat_d, point_num, utmp, vtmp, max_count);

            cudaDeviceSynchronize();
        }
        if (damp_w_to_mean) {
            cudaMemset(wtmp, 0, sizeof(double) * nlat_bins * nv * max_count);
            zonal_w<<<NB, NTH>>>(W_d, Rho_d, zonal_mean_tab_d, point_num, wtmp, max_count);

            cudaDeviceSynchronize();
        }

        if (damp_uv_to_mean && damp_w_to_mean) {
            gpu_sum_on_device_sponge<1024>(utmp, max_count, vbar_h, nv, nlat_bins, 0);
            gpu_sum_on_device_sponge<1024>(vtmp, max_count, vbar_h, nv, nlat_bins, 1);
            gpu_sum_on_device_sponge<1024>(wtmp, max_count, vbar_h, nv, nlat_bins, 2);
        }
        else if (damp_uv_to_mean) {
            gpu_sum_on_device_sponge<1024>(utmp, max_count, vbar_h, nv, nlat_bins, 0);
            gpu_sum_on_device_sponge<1024>(vtmp, max_count, vbar_h, nv, nlat_bins, 1);
        }
        else if (damp_w_to_mean) {
            gpu_sum_on_device_sponge<1024>(wtmp, max_count, vbar_h, nv, nlat_bins, 2);
        }
        if (damp_uv_to_mean || damp_w_to_mean) {
            cudaMemset(vbar_d, 0, sizeof(double) * 3 * nlat_bins * nv);
            cudaMemcpy(vbar_d, vbar_h, 3 * nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);
        }

        if (sim.RayleighSpongeT) {
            zonal_temp<<<NB, NTH>>>(pressure_d,
                                    Rho_d,
                                    Tbar_d,
                                    zonal_mean_tab_d,
                                    lonlat_d,
                                    point_num,
                                    Ttmp,
                                    Rd_d,
                                    max_count);

            cudaDeviceSynchronize();

            int ilat, lev;
            for (ilat = 0; ilat < nlat_bins; ilat++) {
                for (lev = 0; lev < nv; lev++) {
                    Tbar_h[ilat * nv + lev] = gpu_sum_on_device<1024>(
                        &(Ttmp[ilat * nv * max_count + lev * max_count]), max_count);
                }
            }
            cudaMemcpy(Tbar_d, Tbar_h, nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);
        }

        double Rv_fac = 1;
        if (shrink_sponge == true) {
            if (current_step * timestep >= t_shrink * timestep) {
                double shrink_scale = timestep * 1000;
                Rv_fac = exp(-(current_step * timestep - t_shrink * timestep) / shrink_scale);
            }
        }

        sponge_layer<<<NBRT, NTH>>>(Mh_d,
                                    Rho_d,
                                    W_d,
                                    Wh_d,
                                    pressure_d,
                                    vbar_d,
                                    Tbar_d,
                                    zonal_mean_tab_d,
                                    lonlat_d,
                                    Altitude_d,
                                    Altitudeh_d,
                                    Ruv_sponge,
                                    Rw_sponge,
                                    RT_sponge,
                                    Rv_fac,
                                    ns_ray_sponge,
                                    damp_uv_to_mean,
                                    damp_w_to_mean,
                                    true,
                                    timestep,
                                    Rd_d,
                                    nlat_bins,
                                    point_num,
                                    nv,
                                    sim.RayleighSpongeT,
                                    profx_dMh_d,
                                    profx_dWh_d,
                                    profx_dW_d,
                                    profx_Qheat_d);

        BENCH_POINT_I(
            current_step, "phy_Sponge", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "W_d"))
    }


    //  Computes the initial temperature.
    bool      calcT = true;
    if (ultrahot_thermo == VARY_R_CP) {
        check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        update_temperature_Rd_Cp<<<NB, NTH>>>(temperature_d,
                                              Rd_d,
                                              Cp_d,
                                              pressure_d,
                                              Rho_d,
                                              GibbsT_d,
                                              GibbsdG_d,
                                              GibbsN,
                                              point_num,
                                              check_d);

        cudaDeviceSynchronize();
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
        if (check_h) {
            log::printf("\n\n Ridder's method failed for T and Rd\n");
            exit(EXIT_FAILURE);
        }
        calcT = false;
    }
    Compute_temperature<<<NB, NTH>>>(
        temperature_d, pt_d, pressure_d, Rho_d, sim.P_Ref, Rd_d, Cp_d, point_num, calcT);

    // cudaMemcpy(W_h, W_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
    // for (int lev = 0; lev < nv; lev++) {
    //     printf("%d %.15e\n", lev, W_h[0 * nv + lev]);
    // }
    BENCH_POINT_I(
        current_step, "phy_T", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

#ifdef BENCH_NAN_CHECK
    check_h = check_array_for_nan(temperature_d, nv * point_num, 1, check_d);
    if (check_h) {
        log::printf("\n\n Error in NAN check after PROFX:compute_temp!\n");
        exit(EXIT_FAILURE);
    }
#endif

    if (sim.conv_adj) {
        if (current_step > 100) {

            cudaDeviceSynchronize();

            bool ray_mode = true;

            if (ray_mode) {

                printf("start ray_dry_conv_adj \n");

                double *Pressure_d,    // Pressure [Pa]
                             double *Pressureh_d,   // Mid-point pressure [Pa]
                             double *dT_conv_d,
                             double *Temperature_d, // Temperature [K]
                             double *profx_Qheat_d,
                             double *pt_d,          // Potential temperature [K]
                             double *Rho_d,         // Density [m^3/kg]
                             double *Cp_d,          // Specific heat capacity [J/kg/K]
                             double *Rd_d,          // Gas constant [J/kg/K]
                             double  Gravit,        // Gravity [m/s^2]
                             double *Altitude_d,    // Altitudes of the layers
                             double *Altitudeh_d,   // Altitudes of the interfaces
                             double timestep,       // time step [s]
                             int conv_adj_iter, // number of iterations of entire algorithm allowed
                             bool soft_adjust,
                             int num,           // Number of columns
                             int nv)

                ray_dry_conv_adj<<<NBRT, NTH>>>(pressure_d,    // Pressure [Pa]
                                                pressureh_d,   // mid-point pressure [Pa]
                                                dT_conv_d,
                                                temperature_d, // Temperature [K]
                                                profx_Qheat_d,
                                                pt_d,          // Pot temperature [K]
                                                Rho_d,         // Density [m^3/kg]
                                                Cp_d,          // Specific heat capacity [J/kg/K]
                                                Rd_d,          // Gas constant [J/kg/K]
                                                sim.Gravit,    // Gravity [m/s^2]
                                                Altitude_d,    // Altitudes of the layers
                                                Altitudeh_d,   // Altitudes of the interfaces
                                                timestep,      // time step [s]
                                                sim.conv_adj_iter,
                                                sim.soft_adjustment,
                                                point_num, // Number of columns
                                                nv);       // number of vertical layers

                 printf("end ray_dry_conv_adj \n");
            } else {   

                dry_conv_adj<<<NBRT, NTH>>>(pressure_d,    // Pressure [Pa]
                                            pressureh_d,   // mid-point pressure [Pa]
                                            temperature_d, // Temperature [K]
                                            profx_Qheat_d,
                                            pt_d,        // Pot temperature [K]
                                            Rho_d,       // Density [m^3/kg]
                                            Cp_d,        // Specific heat capacity [J/kg/K]
                                            Rd_d,        // Gas constant [J/kg/K]
                                            sim.Gravit,  // Gravity [m/s^2]
                                            Altitude_d,  // Altitudes of the layers
                                            Altitudeh_d, // Altitudes of the interfaces
                                            timestep,
                                            sim.conv_adj_iter,
                                            sim.soft_adjustment,
                                            point_num, // Number of columns
                                            nv);       // number of vertical layers
            }
        }
    }

    BENCH_POINT_I(current_step,
                  "dry_conv_adj ",
                  (),
                  ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))


    ///////////////////////
    // HELD SUAREZ TEST  //
    ///////////////////////
    //
    if (core_benchmark == HELD_SUAREZ) {
        cudaDeviceSynchronize();
        held_suarez<<<NB, NTH>>>(Mh_d,
                                 pressure_d,
                                 Rho_d,
                                 temperature_d,
                                 sim.Gravit,
                                 sim.Cp,
                                 sim.Rd,
                                 Altitude_d,
                                 Altitudeh_d,
                                 lonlat_d,
                                 timestep,
                                 point_num);
    }
    else if (core_benchmark == TIDALLY_LOCKED_EARTH) {
        cudaDeviceSynchronize();
        tidalearth<<<NB, NTH>>>(Mh_d,
                                pressure_d,
                                Rho_d,
                                temperature_d,
                                sim.Gravit,
                                sim.Cp,
                                sim.Rd,
                                Altitude_d,
                                Altitudeh_d,
                                lonlat_d,
                                timestep,
                                point_num);
    }
    else if (core_benchmark == SHALLOW_HOT_JUPITER) {
        cudaDeviceSynchronize();
        shallowHJ<<<NB, NTH>>>(Mh_d,
                               pressure_d,
                               Rho_d,
                               temperature_d,
                               sim.Gravit,
                               sim.Cp,
                               sim.Rd,
                               Altitude_d,
                               Altitudeh_d,
                               lonlat_d,
                               timestep,
                               point_num);
    }
    else if (core_benchmark == DEEP_HOT_JUPITER) {
        cudaDeviceSynchronize();
        deepHJ<<<NB, NTH>>>(Mh_d,
                            pressure_d,
                            Rho_d,
                            temperature_d,
                            profx_Qheat_d,
                            sim.Gravit,
                            sim.Cp,
                            sim.Rd,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            timestep,
                            point_num);
    }
    // else if (core_benchmark == ACOUSTIC_TEST) {
    //     if (current_step == 1) {
    //         acoustic_test<<<NB, NTH>>>(pressure_d,
    //                                    Rho_d,
    //                                    temperature_d,
    //                                    sim.Rd,
    //                                    Altitude_d,
    //                                    lonlat_d,
    //                                    sim.Top_altitude,
    //                                    point_num);
    //     }
    // }
    // else if (core_benchmark == GWAVE_TEST) {
    //     if (current_step == 1) {
    //         gwave_test<<<NB, NTH>>>(pressure_d,
    //                                 Rho_d,
    //                                 temperature_d,
    //                                 sim.Rd,
    //                                 sim.Cp,
    //                                 sim.P_Ref,
    //                                 Altitude_d,
    //                                 lonlat_d,
    //                                 sim.Top_altitude,
    //                                 point_num);
    //     }
    // }
    //always do this nan check so the code doesn't keep computing garbage
    check_h = false;
    cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
    isnan_check<<<16, NTH>>>(temperature_d, nv, point_num, check_d);
    cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
    if (check_h) {
        log::printf("\n\n Error in NAN check before profx::phy_modules_phy_loop!\n");
        exit(EXIT_FAILURE);
    }


    BENCH_POINT_I(current_step,
                  "phy_core_benchmark",
                  (),
                  ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    // here is where we call all "external" physics modules
    if (phy_modules_execute) {
        cudaDeviceSynchronize();
        insolation.phy_loop(*this,
                            sim,
                            current_step, // Step number
                            timestep);    // Time-step [s]
        phy_modules_phy_loop(*this,
                             sim,
                             diag,
                             current_step, // Step number
                             timestep);    // Time-step [s]
    }


    BENCH_POINT_I(current_step,
                  "phy_module",
                  (),
                  ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d", "Qheat"))
    //  Computes the new pressures.
    cudaDeviceSynchronize();

    //apply heating here if gcm_off = true
    if (sim.gcm_off == true) {
        apply_heating<<<NB, NTH>>>(
            temperature_d, profx_Qheat_d, Rho_d, Cp_d, Rd_d, timestep, point_num);
        cudaDeviceSynchronize();
        Compute_pressure_density_hydrostatic<<<NBRT, NTH>>>(pressure_d,
                                                            Rho_d,
                                                            temperature_d,
                                                            Tsurface_d,
                                                            Rd_d,
                                                            Altitude_d,
                                                            sim.P_Ref,
                                                            sim.Gravit,
                                                            point_num,
                                                            nv,
                                                            surface);
    }
    else {
        Compute_pressure<<<NB, NTH>>>(pressure_d, temperature_d, Rho_d, Rd_d, point_num);
    }

    cudaDeviceSynchronize();

    //always do this nan check so the code doesn't keep computing garbage
    check_h = false;
    cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
    isnan_check<<<16, NTH>>>(temperature_d, nv, point_num, check_d);
    cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
    if (check_h) {
        log::printf("\n\n Error in NAN check after PROFX:compute_pressure!\n");
        exit(EXIT_FAILURE);
    }

#ifdef BENCHMARKING
    // recompute temperature from pressure and density, to avoid rounding issues when comparing
    Compute_temperature_only<<<NB, NTH>>>(temperature_d, pressure_d, Rho_d, Rd_d, point_num);
#endif // BENCHMARKING

    BENCH_POINT_I(current_step,
                  "phy_END",
                  (),
                  ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"));

    if (sim.out_interm_momentum) {
        cudaMemcpy(Mh_profx_d, Mh_d, point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    }
}

void ESP::globdiag(const SimulationSetup& sim) {
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3 NB((point_num / NTH) + 1, nv, 1);
    dim3 NB1layer((point_num / NTH + 1), 1, 1);

    // calculate quantities we hope to conserve!
    cudaMemset(GlobalE_d, 0, sizeof(double));
    cudaMemset(GlobalMass_d, 0, sizeof(double));
    cudaMemset(GlobalAMx_d, 0, sizeof(double));
    cudaMemset(GlobalAMy_d, 0, sizeof(double));
    cudaMemset(GlobalAMz_d, 0, sizeof(double));
    cudaMemset(GlobalEnt_d, 0, sizeof(double));

    CalcMass<<<NB, NTH>>>(Mass_d,
                          GlobalMass_d,
                          Rho_d,
                          sim.A,
                          Altitudeh_d,
                          lonlat_d,
                          areasT_d,
                          point_num,
                          sim.DeepModel);

    CalcTotEnergy<<<NB, NTH>>>(Etotal_d,
                               GlobalE_d,
                               Mh_d,
                               W_d,
                               Rho_d,
                               temperature_d,
                               sim.Gravit,
                               Cp_d,
                               Rd_d,
                               sim.A,
                               Altitude_d,
                               Altitudeh_d,
                               lonlat_d,
                               areasT_d,
                               func_r_d,
                               point_num,
                               sim.DeepModel);

    CalcAngMom<<<NB, NTH>>>(AngMomx_d,
                            AngMomy_d,
                            AngMomz_d,
                            GlobalAMx_d,
                            GlobalAMy_d,
                            GlobalAMz_d,
                            Mh_d,
                            Rho_d,
                            sim.A,
                            sim.Omega,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            areasT_d,
                            func_r_d,
                            point_num,
                            sim.DeepModel);

    CalcEntropy<<<NB, NTH>>>(Entropy_d,
                             pressure_d,
                             temperature_d,
                             Cp_d,
                             Rd_d,
                             sim.A,
                             sim.P_Ref,
                             Altitude_d,
                             Altitudeh_d,
                             lonlat_d,
                             areasT_d,
                             func_r_d,
                             point_num,
                             sim.DeepModel);

    // run globdiag on device
    // compute globals
    GlobalE_h    = gpu_sum_on_device<1024>(Etotal_d, point_num * nv);
    GlobalMass_h = gpu_sum_on_device<1024>(Mass_d, point_num * nv);
    GlobalAMx_h  = gpu_sum_on_device<1024>(AngMomx_d, point_num * nv);
    GlobalAMy_h  = gpu_sum_on_device<1024>(AngMomy_d, point_num * nv);
    GlobalAMz_h  = gpu_sum_on_device<1024>(AngMomz_d, point_num * nv);
    GlobalEnt_h  = gpu_sum_on_device<1024>(Entropy_d, point_num * nv);

    if (surface) {
        EnergySurface<<<NB1layer, NTH>>>(Esurf_d, Tsurface_d, areasT_d, Csurf, point_num);
        double GlobalEsurf = 0.0;
        GlobalEsurf        = gpu_sum_on_device<1024>(Esurf_d, point_num);
        GlobalE_h += GlobalEsurf;
    }
}

__global__ void update_mean(double* pressure_mean_d,
                            double* pressure_d,
                            double* Rho_mean_d,
                            double* Rho_d,
                            double* Mh_mean_d,
                            double* Mh_d,
                            double* Wh_mean_d,
                            double* Wh_d,
                            int     n_since_out,
                            int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        pressure_mean_d[id * nv + lev] =
            1.0 / n_since_out
            * (pressure_mean_d[id * nv + lev] * (n_since_out - 1) + pressure_d[id * nv + lev]);
        Rho_mean_d[id * nv + lev] =
            1.0 / n_since_out
            * (Rho_mean_d[id * nv + lev] * (n_since_out - 1) + Rho_d[id * nv + lev]);
        Mh_mean_d[3 * id * nv + 3 * lev + 0] =
            1.0 / n_since_out
            * (Mh_mean_d[3 * id * nv + 3 * lev + 0] * (n_since_out - 1)
               + Mh_d[3 * id * nv + 3 * lev] + 0);
        Mh_mean_d[3 * id * nv + 3 * lev + 1] =
            1.0 / n_since_out
            * (Mh_mean_d[3 * id * nv + 3 * lev + 1] * (n_since_out - 1)
               + Mh_d[3 * id * nv + 3 * lev + 1]);
        Mh_mean_d[3 * id * nv + 3 * lev + 2] =
            1.0 / n_since_out
            * (Mh_mean_d[3 * id * nv + 3 * lev + 2] * (n_since_out - 1)
               + Mh_d[3 * id * nv + 3 * lev + 2]);
        Wh_mean_d[id * (nv + 1) + lev] =
            1.0 / n_since_out
            * (Wh_mean_d[id * (nv + 1) + lev] * (n_since_out - 1) + Wh_d[id * (nv + 1) + lev]);
        if (lev == nv - 1) {
            Wh_mean_d[id * (nv + 1) + lev + 1] =
                1.0 / n_since_out
                * (Wh_mean_d[id * (nv + 1) + lev + 1] * (n_since_out - 1)
                   + Wh_d[id * (nv + 1) + lev + 1]);
        }
    }
}


__host__ void ESP::update_mean_outputs(int n_since_out) {
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3 NB((point_num / NTH) + 1, nv, 1);

    update_mean<<<NB, NTH>>>(pressure_mean_d,
                             pressure_d,
                             Rho_mean_d,
                             Rho_d,
                             Mh_mean_d,
                             Mh_d,
                             Wh_mean_d,
                             Wh_d,
                             n_since_out,
                             point_num);
}
