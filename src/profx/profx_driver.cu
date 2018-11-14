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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "esp.h"
#include "phy/profx_sponge.h"
#include "phy/chemistry_device.h" // Simple chemistry.
#include "phy/dry_conv_adj.h"
#include "phy/profx_auxiliary.h"
#include "phy/profx_deepHJ_hs.h"
#include "phy/profx_held_suarez.h"
#include "phy/profx_shallowHJ_hs.h"
#include "phy/profx_tidalearth_hs.h"
#include "phy/profx_conservation.h"

#include "binary_test.h"
#include "debug_helpers.h"

#include "phy_modules.h"

#include "reduction_add.h"

__host__ void ESP::ProfX(int    chemistry, // Use chemistry
                         int    conv,      //
                         double Omega,     // Rotation rate [1/s]
                         double Cp,        // Specific heat capacity [J/kg/K]
                         double Rd,        // Gas constant [J/kg/K]
                         double mu,        // Atomic mass unit [kg]
                         double kb,        // Boltzmann constant [J/K]
                         double P_Ref,     // Reference pressure [Pa]
                         double Gravit,    // Gravity [m/s^2]
                         double A,         // Planet radius [m]
                         bool   DeepModel,
                         int    n_out,         // output step (triggers conservation calc)
                         bool   sponge,        // Use sponge layer?
                         bool   shrink_sponge, // Shrink sponge after some time (Bonjour Urs!)
                         bool   conservation) {  // calc/output conservation quantities
    USE_BENCHMARK()
    //
    //  Number of threads per block.
    const int NTH = 256;

    //  Specify the block sizes.
    dim3      NB((point_num / NTH) + 1, nv, 1);
    dim3      NBRT((point_num / NTH) + 1, 1, 1);
    dim3      NBTR((point_num / NTH) + 1, nv, ntr);

    if (sponge == true) {
        dim3 NBT((point_num / NTH) + 1, nv, 1);

        cudaMemset(vbar_d, 0, sizeof(double) * 3 * nlat * nv);
        zonal_v<<<NBT, NTH>>>(Mh_d,
                              W_d,
                              Rho_d,
                              vbar_d,
                              zonal_mean_tab_d,
                              lonlat_d,
                              point_num);

        cudaDeviceSynchronize();

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
    BENCH_POINT_I(current_step, "phy_Sponge", vector<string>({}), vector<string>({"Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"}))

    //  Computes the initial temperature.
    Compute_temperature<<<NB, NTH>>>(temperature_d,
                                     pt_d,
                                     pressure_d,
                                     Rho_d,
                                     P_Ref,
                                     Rd,
                                     Cp,
                                     point_num);

    BENCH_POINT_I(current_step, "phy_T", vector<string>({}), vector<string>({"Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"}))

#ifdef BENCH_NAN_CHECK
    check_h = check_array_for_nan(temperature_d, nv * point_num, 1, check_d);
    if (check_h) {
        printf("\n\n Error in NAN check after PROFX:compute_temp!\n");
        exit(EXIT_FAILURE);
    }
#endif

    if (conv) {
        cudaDeviceSynchronize();
        dry_conv_adj<<<NB, NTH>>>(pressure_d,    // Pressure [Pa]
                                  pressureh_d,   // mid-point pressure [Pa]
                                  temperature_d, // Temperature [K]
                                  pt_d,          // Pot temperature [K]
                                  Rho_d,         // Density [m^3/kg]
                                  Cp,            // Specific heat capacity [J/kg/K]
                                  Rd,            // Gas constant [J/kg/K]
                                  Gravit,        // Gravity [m/s^2]
                                  Altitude_d,    // Altitudes of the layers
                                  Altitudeh_d,   // Altitudes of the interfaces
                                  point_num,     // Number of columns
                                  nv);           // number of vertical layers
    }

    BENCH_POINT_I(current_step, "dry_conv_adj ", vector<string>({}), vector<string>({"Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"}))


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
                                 Gravit,
                                 Cp,
                                 Rd,
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
                                Gravit,
                                Cp,
                                Rd,
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
                               Gravit,
                               Cp,
                               Rd,
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
                            Gravit,
                            Cp,
                            Rd,
                            Altitude_d,
                            Altitudeh_d,
                            lonlat_d,
                            timestep,
                            point_num);
    }

    //
    ////////////////////////
    // Simple chemistry
    if (chemistry == 1) {
        cudaDeviceSynchronize();
        Tracers_relax_chemistry_co2<<<NBTR, NTH>>>(tracer_d,
                                                   tauch4_d,
                                                   tauco_d,
                                                   tauh2o_d,
                                                   tauco2_d,
                                                   taunh3_d,
                                                   ch4eq_d,
                                                   coeq_d,
                                                   h2oeq_d,
                                                   co2eq_d,
                                                   nh3eq_d,
                                                   P_che_d,
                                                   T_che_d,
                                                   temperature_d,
                                                   pressure_d,
                                                   Rho_d,
                                                   timestep,
                                                   ntr,
                                                   point_num);
        cudaDeviceSynchronize();
        Tracers_relax_chemistry<<<NBTR, NTH>>>(tracer_d,
                                               tauch4_d,
                                               tauco_d,
                                               tauh2o_d,
                                               tauco2_d,
                                               taunh3_d,
                                               ch4eq_d,
                                               coeq_d,
                                               h2oeq_d,
                                               co2eq_d,
                                               nh3eq_d,
                                               P_che_d,
                                               T_che_d,
                                               temperature_d,
                                               pressure_d,
                                               Rho_d,
                                               timestep,
                                               ntr,
                                               point_num);
    }


    if (core_benchmark == NO_BENCHMARK) {
        cudaDeviceSynchronize();
        phy_modules_mainloop(*this,
                             current_step,   // Step number
                             core_benchmark, // Held-Suarez test option
                             timestep,       // Time-step [s]
                             Omega,          // Rotation rate [1/s]
                             Cp,             // Specific heat capacity [J/kg/K]
                             Rd,             // Gas constant [J/kg/K]
                             mu,             // Atomic mass unit [kg]
                             kb,             // Boltzmann constant [J/K]
                             P_Ref,          // Reference pressure [Pa]
                             Gravit,         // Gravity [m/s^2]
                             A               // Planet radius [m]
        );
    }

    BENCH_POINT_I(current_step, "phy_core_benchmark ", vector<string>({}), vector<string>({"Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"}))
    //  Computes the new pressures.
    cudaDeviceSynchronize();
    Compute_pressure<<<NB, NTH>>>(pressure_d,
                                  temperature_d,
                                  Rho_d,
                                  Rd,
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
                                          Rd,
                                          point_num);
#endif // BENCHMARKING

    BENCH_POINT_I(current_step, "phy_END", vector<string>({}), vector<string>({"Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"}));


    //
    //END OF INTEGRATION
    //
}

// TODO: get constants out of arguments
void ESP::conservation(int    chemistry, //
                       double Omega,     // Rotation rate [1/s]
                       double Cp,        // Specific heat capacity [J/kg/K]
                       double Rd,        // Gas constant [J/kg/K]
                       double mu,        // Atomic mass unit [kg]
                       double kb,        // Boltzmann constant [J/K]
                       double P_Ref,     // Reference pressure [Pa]
                       double Gravit,    // Gravity [m/s^2]
                       double A,         // Planet radius [m]
                       bool   DeepModel) {
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
                          A,
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
                               Gravit,
                               Cp,
                               Rd,
                               A,
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
                            A,
                            Omega,
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
