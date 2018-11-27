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
#include "phy/profx_auxiliary.h"
#include "phy/profx_conservation.h"
#include "phy/profx_deepHJ.h"
#include "phy/profx_held_suarez.h"
#include "phy/profx_shallowHJ.h"
#include "phy/profx_sponge.h"
#include "phy/profx_tidalearth.h"

#include "binary_test.h"
#include "debug_helpers.h"

#include "phy_modules.h"

#include "reduction_add.h"

__host__ void ESP::ProfX(const XPlanet& Planet,
                         int            conv, //
                         bool           DeepModel,
                         int            n_out,         // output step (triggers conservation calc)
                         bool           sponge,        // Use sponge layer?
                         bool           shrink_sponge, // Shrink sponge after some time (Bonjour Urs!)
                         bool           conservation) {          // calc/output conservation quantities
    USE_BENCHMARK()
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3      NB((point_num / NTH) + 1, nv, 1);
    dim3      NBRT((point_num / NTH) + 1, 1, 1);

    if (sponge == true) {
        //dim3 NBT((point_num / NTH) + 1, nv, 1);

        cudaMemset(vbar_d, 0, sizeof(double) * 3 * nlat * nv);
        cudaMemset(utmp, 0, sizeof(double) * nlat * nv * max_count);
        cudaMemset(vtmp, 0, sizeof(double) * nlat * nv * max_count);
        cudaMemset(wtmp, 0, sizeof(double) * nlat * nv * max_count);

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

#ifdef GLOBAL_CONSERVATION_REDUCTIONADD
        int ilat, lev, ind;
        for (ilat = 0; ilat < nlat; ilat++) {
            // vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 0] = 0;
            // vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 1] = 0;
            // vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 2] = 0;
            // for (ind = 0; ind < max_count; ind++) {
            //     vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 0] += utmp[ilat * nv + 3 + (nv - 1) * 3 + ind];
            //     vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 1] += vtmp[ilat * nv + 3 + (nv - 1) * 3 + ind];
            //     vbar_h[ilat * nv + 3 + (nv - 1) * 3 + 2] += wtmp[ilat * nv + 3 + (nv - 1) * 3 + ind];
            // }

            for (lev = 0; lev < nv; lev++) {
                vbar_h[ilat * nv * 3 + lev * 3 + 0] = gpu_sum_on_device<1024>(&(utmp[ilat * nv * max_count + lev * max_count]), max_count);
                vbar_h[ilat * nv * 3 + lev * 3 + 1] = gpu_sum_on_device<1024>(&(vtmp[ilat * nv * max_count + lev * max_count]), max_count);
                vbar_h[ilat * nv * 3 + lev * 3 + 2] = gpu_sum_on_device<1024>(&(wtmp[ilat * nv * max_count + lev * max_count]), max_count);
            }
        }
        cudaMemcpy(vbar_d, vbar_h, 3 * nlat * nv * sizeof(double), cudaMemcpyHostToDevice);
#endif

        // print_vbar(vbar_h, nlat, nv);

        if (shrink_sponge == true) {
            if (current_step * timestep >= t_shrink * 86400) {
                ns_sponge     = 1.0 - 0.5 * (1.0 - ns_sponge);
                shrink_sponge = false;
            }
        }

        sponge_layer<<<NB, NTH>>>(Mh_d,
                                  Rho_d,
                                  W_d,
                                  Wh_d,
                                  vbar_d,
                                  zonal_mean_tab_d,
                                  lonlat_d,
                                  Altitude_d,
                                  Altitudeh_d,
                                  Rv_sponge,
                                  ns_sponge,
                                  timestep,
                                  nlat,
                                  point_num,
                                  nv);
    }
    BENCH_POINT_I(current_step, "phy_Sponge", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    //  Computes the initial temperature.
    Compute_temperature<<<NB, NTH>>>(temperature_d,
                                     pt_d,
                                     pressure_d,
                                     Rho_d,
                                     Planet.P_Ref,
                                     Planet.Rd,
                                     Planet.Cp,
                                     point_num);

    BENCH_POINT_I(current_step, "phy_T", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

#ifdef BENCH_NAN_CHECK
    check_h = check_array_for_nan(temperature_d, nv * point_num, 1, check_d);
    if (check_h) {
        printf("\n\n Error in NAN check after PROFX:compute_temp!\n");
        exit(EXIT_FAILURE);
    }
#endif

    if (conv) {
        cudaDeviceSynchronize();
        dry_conv_adj<<<NBRT, NTH>>>(pressure_d,    // Pressure [Pa]
                                    pressureh_d,   // mid-point pressure [Pa]
                                    temperature_d, // Temperature [K]
                                    pt_d,          // Pot temperature [K]
                                    Rho_d,         // Density [m^3/kg]
                                    Planet.Cp,     // Specific heat capacity [J/kg/K]
                                    Planet.Rd,     // Gas constant [J/kg/K]
                                    Planet.Gravit, // Gravity [m/s^2]
                                    Altitude_d,    // Altitudes of the layers
                                    Altitudeh_d,   // Altitudes of the interfaces
                                    point_num,     // Number of columns
                                    nv);           // number of vertical layers
    }

    BENCH_POINT_I(current_step, "dry_conv_adj ", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))


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
                                 Planet.Gravit,
                                 Planet.Cp,
                                 Planet.Rd,
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
                                Planet.Gravit,
                                Planet.Cp,
                                Planet.Rd,
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
                               Planet.Gravit,
                               Planet.Cp,
                               Planet.Rd,
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
                            Planet.Gravit,
                            Planet.Cp,
                            Planet.Rd,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            timestep,
                            point_num);
    }

    BENCH_POINT_I_PHY(current_step, "phy_core_benchmark", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    if (phy_modules_execute) {
        cudaDeviceSynchronize();
        phy_modules_phy_loop(*this,
                             Planet,
                             current_step, // Step number
                             timestep);    // Time-step [s]
    }

    BENCH_POINT_I_PHY(current_step, "phy_module", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))
    //  Computes the new pressures.
    cudaDeviceSynchronize();
    Compute_pressure<<<NB, NTH>>>(pressure_d,
                                  temperature_d,
                                  Rho_d,
                                  Planet.Rd,
                                  point_num);

    //always do this nan check so the code doesn't keep computing garbage
    check_h = false;
    cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
    isnan_check<<<16, NTH>>>(temperature_d, nv, point_num, check_d);
    cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
    if (check_h) {
        printf("\n\n Error in NAN check after PROFX:compute_pressure!\n");
        exit(EXIT_FAILURE);
    }

#ifdef BENCHMARKING
    // recompute temperature from pressure and density, to avoid rounding issues when comparing
    Compute_temperature_only<<<NB, NTH>>>(temperature_d,
                                          pressure_d,
                                          Rho_d,
                                          Planet.Rd,
                                          point_num);
#endif // BENCHMARKING

    BENCH_POINT_I(current_step, "phy_END", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"));
    
    //
    //END OF INTEGRATION
    //
}

// TODO: get constants out of arguments
void ESP::conservation(const XPlanet& Planet, // planet
                       bool           DeepModel) {
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

    CalcMass<<<NB, NTH>>>(Mass_d,
                          GlobalMass_d,
                          Rho_d,
                          planet.A,
                          Altitudeh_d,
                          lonlat_d,
                          areasT_d,
                          point_num,
                          DeepModel);

    CalcTotEnergy<<<NB, NTH>>>(Etotal_d,
                               GlobalE_d,
                               Mh_d,
                               W_d,
                               Rho_d,
                               temperature_d,
                               Planet.Gravit,
                               Planet.Cp,
                               Planet.Rd,
                               Planet.A,
                               Altitude_d,
                               Altitudeh_d,
                               lonlat_d,
                               areasT_d,
                               func_r_d,
                               point_num,
                               DeepModel);

    CalcAngMom<<<NB, NTH>>>(AngMomx_d,
                            AngMomy_d,
                            AngMomz_d,
                            GlobalAMx_d,
                            GlobalAMy_d,
                            GlobalAMz_d,
                            Mh_d,
                            Rho_d,
                            Planet.A,
                            Planet.Omega,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            areasT_d,
                            point_num,
                            DeepModel);
#ifdef GLOBAL_CONSERVATION_ATOMICADD
    // copy global conservation data to host for output
    copy_global_to_host();
#endif // GLOBAL_CONSERVATION_ATOMICADD

#ifdef GLOBAL_CONSERVATION_REDUCTIONADD
    // run conservation on device

    // compute globals
    GlobalE_h    = gpu_sum_on_device<1024>(Etotal_d, point_num * nv);
    GlobalMass_h = gpu_sum_on_device<1024>(Mass_d, point_num * nv);
    GlobalAMx_h  = gpu_sum_on_device<1024>(AngMomx_d, point_num * nv);
    GlobalAMy_h  = gpu_sum_on_device<1024>(AngMomy_d, point_num * nv);
    GlobalAMz_h  = gpu_sum_on_device<1024>(AngMomz_d, point_num * nv);


#endif // GLOBAL_CONSERVATION_REDUCTIONADD

#ifdef GLOBAL_CONSERVATION_CPUADD

    // copy conservation data to host
    copy_conservation_to_host();

    // compute globals
    GlobalE_h    = cpu_sum<1024>(Etotal_h, point_num * nv);
    GlobalMass_h = cpu_sum<1024>(Mass_h, point_num * nv);
    GlobalAMx_h  = cpu_sum<1024>(AngMomx_h, point_num * nv);
    GlobalAMy_h  = cpu_sum<1024>(AngMomy_h, point_num * nv);
    GlobalAMz_h  = cpu_sum<1024>(AngMomz_h, point_num * nv);

#endif // GLOBAL_CONSERVATION_CPUADD
}
