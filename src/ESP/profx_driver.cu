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
#include "phy/profx_conservation.h"
#include "phy/profx_deepHJ.h"
#include "phy/profx_gwave_test.h"
#include "phy/profx_held_suarez.h"
#include "phy/profx_shallowHJ.h"
#include "phy/profx_sponge.h"
#include "phy/profx_tidalearth.h"

#include "binary_test.h"
#include "debug_helpers.h"
#include "log_writer.h"
#include "phy_modules.h"

#include "reduction_add.h"

__host__ void ESP::ProfX(const SimulationSetup& sim,
                         int                    n_out, // output step (triggers conservation calc)
                         bool                   shrink_sponge) {         // Shrink sponge after some time
    USE_BENCHMARK()
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3      NB((point_num / NTH) + 1, nv, 1);
    dim3      NBRT((point_num / NTH) + 1, 1, 1);

    cudaMemset(profx_dP_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(profx_dMh_d, 0, sizeof(double) * 3 * point_num * nv);
    cudaMemset(profx_dWh_d, 0, sizeof(double) * point_num * nvi);
    cudaMemset(profx_dW_d, 0, sizeof(double) * point_num * nv);


    if (sim.RayleighSponge == true) {
        dim3 NBT((point_num / NTH) + 1, nv, 1);

        cudaMemset(vbar_d, 0, sizeof(double) * 3 * nlat_bins * nv);
        cudaMemset(utmp, 0, sizeof(double) * nlat_bins * nv * max_count);
        cudaMemset(vtmp, 0, sizeof(double) * nlat_bins * nv * max_count);
        cudaMemset(wtmp, 0, sizeof(double) * nlat_bins * nv * max_count);

        zonal_v<<<NB, NTH>>>(Mh_d,
                             W_d,
                             Rho_d,
                             vbar_d,
                             zonal_mean_tab_d,
                             lonlat_d,
                             point_num,
                             utmp,
                             vtmp,
                             wtmp,
                             max_count);

        cudaDeviceSynchronize();

        cudaMemcpy(
            utmp_h, utmp, max_count * nlat_bins * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(
            vtmp_h, vtmp, max_count * nlat_bins * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(
            wtmp_h, wtmp, max_count * nlat_bins * nv * sizeof(double), cudaMemcpyDeviceToHost);

        int ilat, lev;
        for (ilat = 0; ilat < nlat_bins; ilat++) {
            for (lev = 0; lev < nv; lev++) {
                vbar_h[ilat * nv * 3 + lev * 3 + 0] = gpu_sum_on_device<1024>(
                    &(utmp[ilat * nv * max_count + lev * max_count]), max_count);
                vbar_h[ilat * nv * 3 + lev * 3 + 1] = gpu_sum_on_device<1024>(
                    &(vtmp[ilat * nv * max_count + lev * max_count]), max_count);
                vbar_h[ilat * nv * 3 + lev * 3 + 2] = gpu_sum_on_device<1024>(
                    &(wtmp[ilat * nv * max_count + lev * max_count]), max_count);
            }
        }
        cudaMemcpy(vbar_d, vbar_h, 3 * nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);

        if (sim.RayleighSpongeT) {
            zonal_temp<<<NB, NTH>>>(pressure_d,
                                    Rho_d,
                                    Tbar_d,
                                    zonal_mean_tab_d,
                                    lonlat_d,
                                    point_num,
                                    Ttmp,
                                    sim.Rd,
                                    max_count);

            cudaDeviceSynchronize();

            cudaMemcpy(
                Ttmp_h, Ttmp, max_count * nlat_bins * nv * sizeof(double), cudaMemcpyDeviceToHost);

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
                // printf("%d, %f\n", current_step, Rv_fac);
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
                                    timestep,
                                    Rd_d,
                                    nlat_bins,
                                    point_num,
                                    nv,
                                    sim.RayleighSpongeT,
                                    profx_dMh_d,
                                    profx_dWh_d,
                                    profx_dW_d);
    }
    BENCH_POINT_I(current_step,
                  "phy_Sponge",
                  (),
                  ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    //  Computes the initial temperature.
    Compute_temperature<<<NB, NTH>>>(
        temperature_d, pt_d, pressure_d, Rho_d, sim.P_Ref, Rd_d, Cp_d, point_num);

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
        cudaDeviceSynchronize();

        dry_conv_adj<<<NBRT, NTH>>>(pressure_d,    // Pressure [Pa]
                                    pressureh_d,   // mid-point pressure [Pa]
                                    temperature_d, // Temperature [K]
                                    pt_d,          // Pot temperature [K]
                                    Rho_d,         // Density [m^3/kg]
                                    Cp_d,          // Specific heat capacity [J/kg/K]
                                    Rd_d,          // Gas constant [J/kg/K]
                                    sim.Gravit,    // Gravity [m/s^2]
                                    Altitude_d,    // Altitudes of the layers
                                    Altitudeh_d,   // Altitudes of the interfaces
                                    point_num,     // Number of columns
                                    nv);           // number of vertical layers
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
                            profx_dP_d,
                            sim.Gravit,
                            sim.Cp,
                            sim.Rd,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            timestep,
                            point_num);
    }
    else if (core_benchmark == ACOUSTIC_TEST) {
        if (current_step == 1) {
            acoustic_test<<<NB, NTH>>>(pressure_d,
                                       Rho_d,
                                       temperature_d,
                                       sim.Rd,
                                       Altitude_d,
                                       lonlat_d,
                                       sim.Top_altitude,
                                       point_num);
        }
    }
    else if (core_benchmark == GWAVE_TEST) {
        if (current_step == 1) {
            gwave_test<<<NB, NTH>>>(pressure_d,
                                    Rho_d,
                                    temperature_d,
                                    sim.Rd,
                                    sim.Cp,
                                    sim.P_Ref,
                                    Altitude_d,
                                    lonlat_d,
                                    sim.Top_altitude,
                                    point_num);
        }
    }
    BENCH_POINT_I_PHY(current_step,
                      "phy_core_benchmark",
                      (),
                      ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    // here is where we call all "external" physics modules
    if (phy_modules_execute) {
        cudaDeviceSynchronize();
        phy_modules_phy_loop(*this,
                             sim,
                             current_step, // Step number
                             timestep);    // Time-step [s]
    }

    BENCH_POINT_I_PHY(current_step,
                      "phy_module",
                      (),
                      ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))
    //  Computes the new pressures.
    cudaDeviceSynchronize();
    Compute_pressure<<<NB, NTH>>>(pressure_d, temperature_d, Rho_d, Rd_d, point_num);

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

    //
    //END OF INTEGRATION
    //
}

void ESP::conservation(const SimulationSetup& sim) {
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3 NB((point_num / NTH) + 1, nv, 1);

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

    // run conservation on device
    // compute globals
    GlobalE_h    = gpu_sum_on_device<1024>(Etotal_d, point_num * nv);
    GlobalMass_h = gpu_sum_on_device<1024>(Mass_d, point_num * nv);
    GlobalAMx_h  = gpu_sum_on_device<1024>(AngMomx_d, point_num * nv);
    GlobalAMy_h  = gpu_sum_on_device<1024>(AngMomy_d, point_num * nv);
    GlobalAMz_h  = gpu_sum_on_device<1024>(AngMomz_d, point_num * nv);
    GlobalEnt_h  = gpu_sum_on_device<1024>(Entropy_d, point_num * nv);
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
