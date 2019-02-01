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
// Description: DYNAMICAL CORE INTEGRATION
//
// Method: RK3 Method & Forward-Backward Method
//
// Known limitations:
//   - It does not include a shock capture scheme.
//
// Known issues:
//   - Operational in just one GPU.
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

#include "esp.h" // Global parameters.

#include "binary_test.h"
#include "debug_helpers.h"
#include "dyn/thor_adv_cor.h" // Advection term.
#include "dyn/thor_auxiliary.h" // Temperature, interal energy, potential tempareture and effective gravity.
#include "dyn/thor_diff.h"         // Hyper-diffusion.
#include "dyn/thor_div.h"          // Divergence damping.
#include "dyn/thor_fastmodes.h"    // Fast terms.
#include "dyn/thor_slowmodes.h"    // Slow terms.
#include "dyn/thor_vertical_int.h" // Vertical momentum.
#include "log_writer.h"

#include "phy_modules.h"

__host__ void ESP::Thor(const SimulationSetup& sim) {
    const int NTH = 256;

    // Vertical Eq only works on vertical stack of data, can run independently, only uses shared
    // memory for intermediate data that is not shared with neighbours.
    // Need to set the block size so that the internal arrays fit in shared memory for each block
    const int num_th_vertical_eq = 32;
    //  Specify the block sizes.
    const int LN = 16;                     // Size of the inner region side.
    dim3      NT(nl_region, nl_region, 1); // Number of threads in a block.
    dim3      NB(nr, nv, 1);               // Number of blocks.
    dim3      NBD(nr, nv, 6);              // Number of blocks in the diffusion routine.
    dim3      NBDP(2, nv, 6);              // Number of blocks in the diffusion routine. (POLES)
    dim3      NBP(2, nv, 1);               // Number of blocks. (POLES)

    dim3 NBALL((point_num / NTH) + 1, nv, 1); //Number of blocks to execute on all grid points


    //  Number of Small steps
    double ns_totald = 6; // Maximum number of small steps in a large step (double ).
    int    ns_totali = 6; // Maximum number of small steps in a large step (integer).
    int    ns_it;         // Number of small steps in each large step.
    double times;         // Sub-timestep.

    //  Initialize local variables used for the time integration.
    cudaDeviceSynchronize();
    cudaMemcpy(Mhk_d, Mh_d, point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Whk_d, Wh_d, point_num * nvi * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Wk_d, W_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Rhok_d, Rho_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(pressurek_d, pressure_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);

    cudaMemset(Mhs_d, 0, sizeof(double) * 3 * point_num * nv);
    cudaMemset(Rhos_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(Whs_d, 0, sizeof(double) * point_num * nvi);
    cudaMemset(Ws_d, 0, sizeof(double) * point_num * nv);
    cudaMemset(pressures_d, 0, sizeof(double) * point_num * nv);

    if (phy_modules_execute) phy_modules_dyn_core_loop_init(*this);

    USE_BENCHMARK();
    BENCH_POINT_I(current_step,
                  "thor_init",
                  (),
                  ("Rho_d",
                   "pressure_d",
                   "Mh_d",
                   "Wh_d",
                   "temperature_d",
                   "W_d" /*, "tracer_d", "tracers_d", "tracerk_d"*/));


    //  Loop for large time integration.
    for (int rk = 0; rk < 3; rk++) {
        //      Local variables to define the length (times) and the number of the small steps (ns_it).
        if (rk == 0) ns_it = 1;
        if (rk == 1) ns_it = ns_totali / 2;
        if (rk == 2) ns_it = ns_totali;

        if (rk == 0) times = timestep / 3.0;
        if (rk == 1) times = timestep / ns_totald;
        if (rk == 2) times = timestep / ns_totald;

        // initialise some memory

        //
        //      Compute advection and coriolis terms.
        cudaMemset(Adv_d, 0, sizeof(double) * 3 * point_num * nv); // Sets every value of Adv_d to
                                                                   // zero.
        cudaDeviceSynchronize();


        // Updates: Adv_d, v_d

        Compute_Advec_Cori1<LN, LN><<<NB, NT>>>((double3*)Adv_d,
                                                (double3*)v_d,
                                                (double3*)Mhk_d,
                                                (double3*)div_d,
                                                Wk_d,
                                                Rhok_d,
                                                Altitude_d,
                                                sim.A,
                                                (double3*)func_r_d,
                                                maps_d,
                                                nl_region,
                                                sim.DeepModel);
        // Updates: Adv_d, v_d
        Compute_Advec_Cori_Poles<6><<<2, 1>>>(Adv_d,
                                              v_d,
                                              Mhk_d,
                                              div_d,
                                              Wk_d,
                                              Rhok_d,
                                              Altitude_d,
                                              sim.A,
                                              func_r_d,
                                              point_local_d,
                                              point_num,
                                              nv,
                                              sim.DeepModel);

        cudaDeviceSynchronize();
        // Updates: Adv_d
        Compute_Advec_Cori2<<<(point_num / NTH) + 1, NTH>>>(Adv_d,
                                                            v_d,
                                                            Whk_d,
                                                            Rhok_d,
                                                            Altitude_d,
                                                            Altitudeh_d,
                                                            sim.Omega,
                                                            sim.A,
                                                            nv,
                                                            point_num,
                                                            sim.DeepModel);

        //
        //      Computes temperature, internal energy, potential temperature and effective gravity.
        cudaDeviceSynchronize();

        BENCH_POINT_I_S(
            current_step,
            rk,
            "Compute_Advec_Cori",
            (),
            ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d", "Adv_d", "v_d"))

        // Updates: temperature_d, h_d, hh_d, pt_d, pth_d, gtil_d, gtilh_d
        Compute_Temperature_H_Pt_Geff<<<(point_num / NTH) + 1, NTH>>>(temperature_d,
                                                                      pressurek_d,
                                                                      Rhok_d,
                                                                      h_d,
                                                                      hh_d,
                                                                      pt_d,
                                                                      pth_d,
                                                                      gtil_d,
                                                                      gtilh_d,
                                                                      Whk_d,
                                                                      sim.P_Ref,
                                                                      sim.Gravit,
                                                                      sim.Cp,
                                                                      sim.Rd,
                                                                      Altitude_d,
                                                                      Altitudeh_d,
                                                                      point_num,
                                                                      nv);


        //      Initializes slow terms.
        cudaDeviceSynchronize();

        BENCH_POINT_I_S(current_step,
                        rk,
                        "Compute_Temperature_H_Pt_Geff",
                        (),
                        ("temperature_d", "h_d", "hh_d", "pt_d", "pth_d", "gtil_d", "gtilh_d"))

        cudaMemset(SlowMh_d, 0, sizeof(double) * 3 * point_num * nv);
        cudaMemset(SlowWh_d, 0, sizeof(double) * point_num * nvi);
        cudaMemset(SlowRho_d, 0, sizeof(double) * point_num * nv);
        cudaMemset(Slowpressure_d, 0, sizeof(double) * point_num * nv);
        //
        //      Hyper-Diffusion.
        if (sim.HyDiff) {
            cudaMemset(diff_d, 0, sizeof(double) * 6 * point_num * nv);
            cudaDeviceSynchronize();
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op<LN, LN><<<NBD, NT>>>(diffmh_d,
                                              diffw_d,
                                              diffrh_d,
                                              diffpr_d,
                                              diff_d,
                                              Mhk_d,
                                              Rhok_d,
                                              temperature_d,
                                              Wk_d,
                                              areasTr_d,
                                              nvecoa_d,
                                              nvecti_d,
                                              nvecte_d,
                                              func_r_d,
                                              Kdh4_d,
                                              Altitude_d,
                                              sim.A,
                                              sim.Rd,
                                              maps_d,
                                              nl_region,
                                              0,
                                              sim.DeepModel);
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op_Poles<5><<<NBDP, 1>>>(diffmh_d,
                                               diffw_d,
                                               diffrh_d,
                                               diffpr_d,
                                               diff_d,
                                               Mhk_d,
                                               Rhok_d,
                                               temperature_d,
                                               Wk_d,
                                               func_r_d,
                                               areasTr_d,
                                               nvecoa_d,
                                               nvecti_d,
                                               nvecte_d,
                                               Kdh4_d,
                                               Altitude_d,
                                               Altitudeh_d,
                                               sim.A,
                                               sim.Rd,
                                               point_local_d,
                                               point_num,
                                               0,
                                               sim.DeepModel);
            cudaDeviceSynchronize();
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op<LN, LN><<<NBD, NT>>>(diffmh_d,
                                              diffw_d,
                                              diffrh_d,
                                              diffpr_d,
                                              diff_d,
                                              Mhk_d,
                                              Rhok_d,
                                              temperature_d,
                                              Wk_d,
                                              areasTr_d,
                                              nvecoa_d,
                                              nvecti_d,
                                              nvecte_d,
                                              func_r_d,
                                              Kdh4_d,
                                              Altitude_d,
                                              sim.A,
                                              sim.Rd,
                                              maps_d,
                                              nl_region,
                                              1,
                                              sim.DeepModel);
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op_Poles<5><<<NBDP, 1>>>(diffmh_d,
                                               diffw_d,
                                               diffrh_d,
                                               diffpr_d,
                                               diff_d,
                                               Mhk_d,
                                               Rhok_d,
                                               temperature_d,
                                               Wk_d,
                                               func_r_d,
                                               areasTr_d,
                                               nvecoa_d,
                                               nvecti_d,
                                               nvecte_d,
                                               Kdh4_d,
                                               Altitude_d,
                                               Altitudeh_d,
                                               sim.A,
                                               sim.Rd,
                                               point_local_d,
                                               point_num,
                                               1,
                                               sim.DeepModel);

            Diffusion_Op_Vert<<<NBALL, NTH>>>(diffmv_d,
                                              diffwv_d,
                                              diffrv_d,
                                              diffprv_d,
                                              diffv_d,
                                              Mhk_d,
                                              Rhok_d,
                                              temperature_d,
                                              Wk_d,
                                              func_r_d,
                                              Kdv6_d,
                                              Altitude_d,
                                              sim.A,
                                              sim.Rd,
                                              point_num,
                                              0,
                                              sim.DeepModel);

            cudaDeviceSynchronize();
            BENCH_POINT_I_S_PHY(current_step,
                                rk,
                                "Diffusion_Op_Vert1",
                                (),
                                ("diffmh_d",
                                 "diffw_d",
                                 "diffrh_d",
                                 "diffpr_d",
                                 "diff_d",
                                 "diffmv_d",
                                 "diffwv_d",
                                 "diffrv_d",
                                 "diffprv_d",
                                 "diffv_d"))

            Diffusion_Op_Vert<<<NBALL, NTH>>>(diffmv_d,
                                              diffwv_d,
                                              diffrv_d,
                                              diffprv_d,
                                              diffv_d,
                                              Mhk_d,
                                              Rhok_d,
                                              temperature_d,
                                              Wk_d,
                                              func_r_d,
                                              Kdv6_d,
                                              Altitude_d,
                                              sim.A,
                                              sim.Rd,
                                              point_num,
                                              1,
                                              sim.DeepModel);
            cudaDeviceSynchronize();

            BENCH_POINT_I_S_PHY(current_step,
                                rk,
                                "Diffusion_Op_Vert2",
                                (),
                                ("diffmh_d",
                                 "diffw_d",
                                 "diffrh_d",
                                 "diffpr_d",
                                 "diff_d",
                                 "diffmv_d",
                                 "diffwv_d",
                                 "diffrv_d",
                                 "diffprv_d",
                                 "diffv_d"))

            Diffusion_Op_Vert<<<NBALL, NTH>>>(diffmv_d,
                                              diffwv_d,
                                              diffrv_d,
                                              diffprv_d,
                                              diffv_d,
                                              Mhk_d,
                                              Rhok_d,
                                              temperature_d,
                                              Wk_d,
                                              func_r_d,
                                              Kdv6_d,
                                              Altitude_d,
                                              sim.A,
                                              sim.Rd,
                                              point_num,
                                              2,
                                              sim.DeepModel);

            cudaDeviceSynchronize();
            Correct_Horizontal<<<NBALL, NTH>>>(diffmh_d, diffmv_d, func_r_d, point_num);

            BENCH_POINT_I_S_PHY(current_step,
                                rk,
                                "Diffusion_Op_Vert",
                                (),
                                ("diffmh_d",
                                 "diffw_d",
                                 "diffrh_d",
                                 "diffpr_d",
                                 "diff_d",
                                 "diffmv_d",
                                 "diffwv_d",
                                 "diffrv_d",
                                 "diffprv_d",
                                 "diffv_d"))
        }

        if (phy_modules_execute)
            phy_modules_dyn_core_loop_slow_modes(*this, sim, current_step, times);

        BENCH_POINT_I_S_PHY(current_step,
                            rk,
                            "DivDamp",
                            (),
                            ("Rhos_d",
                             "Rhok_d",
                             "Mhs_d",
                             "Mhk_d",
                             "Whs_d",
                             "Whk_d",
                             "pressures_d",
                             "pressurek_d",
                             "pressure_d"))
        //
        //      Divergence damping
        cudaMemset(DivM_d, 0, sizeof(double) * point_num * 3 * nv);
        cudaMemset(divg_Mh_d, 0, sizeof(double) * point_num * 3 * nv);
        if (sim.DivDampP) {
            cudaDeviceSynchronize();
            // Updates: DivM_d, divg_Mh_d
            DivM_Op<LN, LN><<<NB, NT>>>(DivM_d,
                                        divg_Mh_d,
                                        Mhk_d,
                                        Whk_d,
                                        Kdhz_d,
                                        areasTr_d,
                                        nvecoa_d,
                                        nvecti_d,
                                        nvecte_d,
                                        func_r_d,
                                        Altitudeh_d,
                                        Altitude_d,
                                        sim.A,
                                        maps_d,
                                        nl_region,
                                        0,
                                        sim.DeepModel);
            // Updates: DivM_d, divg_Mh_d
            DivM_Op_Poles<5><<<NBP, 1>>>(DivM_d,
                                         divg_Mh_d,
                                         Mhk_d,
                                         Whk_d,
                                         Kdhz_d,
                                         areasTr_d,
                                         nvecoa_d,
                                         nvecti_d,
                                         nvecte_d,
                                         func_r_d,
                                         Altitudeh_d,
                                         Altitude_d,
                                         sim.A,
                                         point_local_d,
                                         point_num,
                                         0,
                                         sim.DeepModel);

            cudaDeviceSynchronize();
            // Updates: DivM_d, divg_Mh_d
            DivM_Op<LN, LN><<<NB, NT>>>(DivM_d,
                                        divg_Mh_d,
                                        Mhk_d,
                                        Whk_d,
                                        Kdhz_d,
                                        areasTr_d,
                                        nvecoa_d,
                                        nvecti_d,
                                        nvecte_d,
                                        func_r_d,
                                        Altitudeh_d,
                                        Altitude_d,
                                        sim.A,
                                        maps_d,
                                        nl_region,
                                        1,
                                        sim.DeepModel);
            // Updates: DivM_d, divg_Mh_d
            DivM_Op_Poles<5><<<NBP, 1>>>(DivM_d,
                                         divg_Mh_d,
                                         Mhk_d,
                                         Whk_d,
                                         Kdhz_d,
                                         areasTr_d,
                                         nvecoa_d,
                                         nvecti_d,
                                         nvecte_d,
                                         func_r_d,
                                         Altitudeh_d,
                                         Altitude_d,
                                         sim.A,
                                         point_local_d,
                                         point_num,
                                         1,
                                         sim.DeepModel);
        }

        BENCH_POINT_I_S(current_step, rk, "DivM_Op_Poles", (), ("DivM_d", "divg_Mh_d"))

        //
        //      Slow Modes
        cudaDeviceSynchronize();


        // Updates: SlowMh_d, SlowWh_d, SlowRho_d, Slowpressure_d
        Compute_Slow_Modes<LN, LN><<<NB, NT>>>(SlowMh_d,
                                               SlowWh_d,
                                               SlowRho_d,
                                               Slowpressure_d,
                                               Mhk_d,
                                               Whk_d,
                                               Rhok_d,
                                               Adv_d,
                                               DivM_d,
                                               diffmh_d,
                                               diffw_d,
                                               diffrh_d,
                                               diffpr_d,
                                               diffmv_d,
                                               diffwv_d,
                                               diffrv_d,
                                               diffprv_d,
                                               pressurek_d,
                                               h_d,
                                               hh_d,
                                               gtil_d,
                                               grad_d,
                                               div_d,
                                               Altitude_d,
                                               Altitudeh_d,
                                               sim.A,
                                               sim.Gravit,
                                               sim.Cp,
                                               sim.Rd,
                                               func_r_d,
                                               maps_d,
                                               nl_region,
                                               sim.DeepModel,
                                               sim.NonHydro);
        cudaDeviceSynchronize();
        // Updates: SlowMh_d, SlowWh_d, SlowRho_d, Slowpressure_d
        Compute_Slow_Modes_Poles<6><<<2, 1>>>(SlowMh_d,
                                              SlowWh_d,
                                              SlowRho_d,
                                              Slowpressure_d,
                                              Mhk_d,
                                              Whk_d,
                                              Rhok_d,
                                              Adv_d,
                                              DivM_d,
                                              diffmh_d,
                                              diffw_d,
                                              diffrh_d,
                                              diffpr_d,
                                              diffmv_d,
                                              diffwv_d,
                                              diffrv_d,
                                              diffprv_d,
                                              pressurek_d,
                                              h_d,
                                              hh_d,
                                              gtil_d,
                                              grad_d,
                                              div_d,
                                              Altitude_d,
                                              Altitudeh_d,
                                              sim.A,
                                              sim.Gravit,
                                              sim.Cp,
                                              sim.Rd,
                                              func_r_d,
                                              point_local_d,
                                              nv,
                                              point_num,
                                              sim.DeepModel,
                                              sim.NonHydro);


        BENCH_POINT_I_S(current_step,
                        rk,
                        "Compute_Slow_Modes_Poles",
                        (),
                        ("SlowMh_d", "SlowWh_d", "SlowRho_d", "Slowpressure_d"))

        //      Updates or initializes deviations.
        if (rk > 0) {
            cudaDeviceSynchronize();

            BENCH_POINT_I_S_PHY(current_step,
                                rk,
                                "bRK",
                                (),
                                ("Rhos_d",
                                 "Rhok_d",
                                 "Mhs_d",
                                 "Mhk_d",
                                 "Whs_d",
                                 "Whk_d",
                                 "pressures_d",
                                 "pressurek_d",
                                 "pressure_d"))

            // Updates: Mhs_d, Whs_d, Ws_d, Rhos_d, pressures_d
            UpdateRK<<<(point_num / NTH) + 1, NTH>>>(Mhs_d,
                                                     Mhk_d,
                                                     Mh_d,
                                                     Whs_d,
                                                     Whk_d,
                                                     Wh_d,
                                                     Ws_d,
                                                     Rhos_d,
                                                     Rhok_d,
                                                     Rho_d,
                                                     pressures_d,
                                                     pressurek_d,
                                                     pressure_d,
                                                     func_r_d,
                                                     Altitude_d,
                                                     Altitudeh_d,
                                                     point_num,
                                                     nv);

            BENCH_POINT_I_S_PHY(current_step,
                                rk,
                                "RK",
                                (),
                                ("Rhos_d",
                                 "Rhok_d",
                                 "Mhs_d",
                                 "Mhk_d",
                                 "Whs_d",
                                 "Whk_d",
                                 "pressures_d",
                                 "pressurek_d",
                                 "pressure_d"))
        }

        //
        //      SMALL-STEPS
        // printf("\n << nlarge = %d <<<<<<<<<<<<\n",rk);
        for (int ns = 0; ns < ns_it; ns++) {
            // printf("// nsmall = %d //////////////////\n",ns);
            cudaMemset(DivM_d, 0, sizeof(double) * point_num * 3 * nv);
            cudaMemset(divg_Mh_d, 0, sizeof(double) * point_num * 3 * nv);
            if (sim.DivDampP) {
                cudaDeviceSynchronize();
                // Updates: DivM_d, divg_Mh_d
                DivM_Op<LN, LN><<<NB, NT>>>(DivM_d,
                                            divg_Mh_d,
                                            Mhs_d,
                                            Whs_d,
                                            Kdhz_d,
                                            areasTr_d,
                                            nvecoa_d,
                                            nvecti_d,
                                            nvecte_d,
                                            func_r_d,
                                            Altitudeh_d,
                                            Altitude_d,
                                            sim.A,
                                            maps_d,
                                            nl_region,
                                            0,
                                            sim.DeepModel);
                // Updates: DivM_d, divg_Mh_d
                DivM_Op_Poles<5><<<NBP, 1>>>(DivM_d,
                                             divg_Mh_d,
                                             Mhs_d,
                                             Whs_d,
                                             Kdhz_d,
                                             areasTr_d,
                                             nvecoa_d,
                                             nvecti_d,
                                             nvecte_d,
                                             func_r_d,
                                             Altitudeh_d,
                                             Altitude_d,
                                             sim.A,
                                             point_local_d,
                                             point_num,
                                             0,
                                             sim.DeepModel);

                cudaDeviceSynchronize();
                // Updates: DivM_d, divg_Mh_d
                DivM_Op<LN, LN><<<NB, NT>>>(DivM_d,
                                            divg_Mh_d,
                                            Mhs_d,
                                            Whs_d,
                                            Kdhz_d,
                                            areasTr_d,
                                            nvecoa_d,
                                            nvecti_d,
                                            nvecte_d,
                                            func_r_d,
                                            Altitudeh_d,
                                            Altitude_d,
                                            sim.A,
                                            maps_d,
                                            nl_region,
                                            1,
                                            sim.DeepModel);
                // Updates: DivM_d, divg_Mh_d
                DivM_Op_Poles<5><<<NBP, 1>>>(DivM_d,
                                             divg_Mh_d,
                                             Mhs_d,
                                             Whs_d,
                                             Kdhz_d,
                                             areasTr_d,
                                             nvecoa_d,
                                             nvecti_d,
                                             nvecte_d,
                                             func_r_d,
                                             Altitudeh_d,
                                             Altitude_d,
                                             sim.A,
                                             point_local_d,
                                             point_num,
                                             1,
                                             sim.DeepModel);

                BENCH_POINT_I_SS(current_step, rk, ns, "DivM_Op_Poles", (), ("DivM_d", "divg_Mh_d"))
            }

            //          Momentum equation.
            cudaDeviceSynchronize();
            // Updates: Mhs_d
            Momentum_Eq<LN, LN><<<NB, NT>>>(Mhs_d,
                                            pressures_d,
                                            SlowMh_d,
                                            grad_d,
                                            Altitude_d,
                                            DivM_d,
                                            sim.A,
                                            func_r_d,
                                            times,
                                            maps_d,
                                            nl_region,
                                            sim.DeepModel);
            // Updates: Mhs_d
            Momentum_Eq_Poles<6><<<2, 1>>>(Mhs_d,
                                           pressures_d,
                                           SlowMh_d,
                                           grad_d,
                                           Altitude_d,
                                           DivM_d,
                                           sim.A,
                                           func_r_d,
                                           times,
                                           point_local_d,
                                           nv,
                                           point_num,
                                           sim.DeepModel);

            //          Vertical Momentum
            cudaDeviceSynchronize();

            BENCH_POINT_I_SS(current_step,
                             rk,
                             ns,
                             "Momentum_Eq",
                             (),
                             ("Rho_d", "pressures_d", "Mhs_d", "Wh_d", "temperature_d", "W_d"))
            // Updates: Sp_d, Sd_d
            Prepare_Implicit_Vertical<LN, LN><<<NB, NT>>>(Mhs_d,
                                                          h_d,
                                                          div_d,
                                                          Slowpressure_d,
                                                          SlowRho_d,
                                                          Sp_d,
                                                          Sd_d,
                                                          Altitude_d,
                                                          sim.Cp,
                                                          sim.Rd,
                                                          sim.A,
                                                          maps_d,
                                                          nl_region,
                                                          sim.DeepModel);


            cudaDeviceSynchronize();
            // Updates: Sp_d, Sd_d
            Prepare_Implicit_Vertical_Poles<6><<<2, 1>>>(Mhs_d,
                                                         h_d,
                                                         div_d,
                                                         Slowpressure_d,
                                                         SlowRho_d,
                                                         Sp_d,
                                                         Sd_d,
                                                         Altitude_d,
                                                         sim.Cp,
                                                         sim.Rd,
                                                         sim.A,
                                                         point_local_d,
                                                         point_num,
                                                         nv,
                                                         sim.DeepModel);

            BENCH_POINT_I_SS(
                current_step, rk, ns, "Prepare_Implicit_Vertical", (), ("Sp_d", "Sd_d"))

            cudaDeviceSynchronize();

            // Updates: Whs_d, Ws_d
            Vertical_Eq<<<(point_num / num_th_vertical_eq) + 1,
                          num_th_vertical_eq,
                          2 * num_th_vertical_eq * nvi * sizeof(double)>>>(Whs_d,
                                                                           Ws_d,
                                                                           pressures_d,
                                                                           h_d,
                                                                           hh_d,
                                                                           Rhos_d,
                                                                           gtil_d,
                                                                           gtilh_d,
                                                                           Sp_d,
                                                                           Sd_d,
                                                                           SlowWh_d,
                                                                           sim.Cp,
                                                                           sim.Rd,
                                                                           times,
                                                                           sim.Gravit,
                                                                           Altitude_d,
                                                                           Altitudeh_d,
                                                                           sim.A,
                                                                           sim.NonHydro,
                                                                           point_num,
                                                                           nv,
                                                                           nvi,
                                                                           sim.DeepModel);

            cudaError_t err = cudaGetLastError();


            // Check device query
            if (err != cudaSuccess) {
                log::printf("[%s:%d] CUDA error check reports error: %s\n",
                            __FILE__,
                            __LINE__,
                            cudaGetErrorString(err));
            }

            //          Pressure and density equations.
            cudaDeviceSynchronize();
            BENCH_POINT_I_SS_PHY(current_step,
                                 rk,
                                 ns,
                                 "Vertical_Eq",
                                 (),
                                 ("Whs_d", "Ws_d", "pressures_d", "h_d", "hh_d", "Rhos_d"));

            // update the physics modules in fast mode
            if (phy_modules_execute)
                phy_modules_dyn_core_loop_fast_modes(*this, sim, current_step, times);


            BENCH_POINT_I_SS_PHY(current_step,
                                 rk,
                                 ns,
                                 "Phy_mod_fast_mode",
                                 (),
                                 ("Whs_d", "Ws_d", "pressures_d", "h_d", "hh_d", "Rhos_d"))

            // Updates: pressures_d, Rhos_d
            Density_Pressure_Eqs<LN, LN><<<NB, NT>>>(pressures_d,
                                                     pressurek_d,
                                                     Rhos_d,
                                                     Rhok_d,
                                                     Mhs_d,
                                                     Mhk_d,
                                                     Whs_d,
                                                     Whk_d,
                                                     pt_d,
                                                     pth_d,
                                                     SlowRho_d,
                                                     diffpr_d,
                                                     diffprv_d,
                                                     div_d,
                                                     Altitude_d,
                                                     Altitudeh_d,
                                                     sim.Cp,
                                                     sim.Rd,
                                                     sim.A,
                                                     sim.P_Ref,
                                                     times,
                                                     maps_d,
                                                     nl_region,
                                                     sim.DeepModel);

            cudaDeviceSynchronize();
            BENCH_POINT_I_SS(
                current_step, rk, ns, "Density_Pressure_Eqs", (), ("pressures_d", "Rhos_d"))
            // Updates: pressures_d, Rhos_d
            Density_Pressure_Eqs_Poles<6><<<2, 1>>>(pressures_d,
                                                    pressurek_d,
                                                    Rhos_d,
                                                    Rhok_d,
                                                    Mhs_d,
                                                    Mhk_d,
                                                    Whs_d,
                                                    Whk_d,
                                                    pt_d,
                                                    pth_d,
                                                    SlowRho_d,
                                                    diffpr_d,
                                                    diffprv_d,
                                                    div_d,
                                                    Altitude_d,
                                                    Altitudeh_d,
                                                    sim.Cp,
                                                    sim.Rd,
                                                    sim.A,
                                                    sim.P_Ref,
                                                    times,
                                                    point_local_d,
                                                    point_num,
                                                    nv,
                                                    sim.DeepModel);

            BENCH_POINT_I_SS(
                current_step, rk, ns, "Density_Pressure_Eqs_Poles", (), ("pressures_d", "Rhos_d"))
        }
        BENCH_POINT_I_S_PHY(
            current_step,
            rk,
            "bRK2",
            (),
            ("Rhos_d", "Rhok_d", "Mhs_d", "Mhk_d", "Whs_d", "Whk_d", "pressures_d", "pressurek_d"))
        //      Update quantities for the long loop.
        cudaDeviceSynchronize();
        // Updates: Mhk_d, Whk_d, Wk_d, Rhok_d, pressurek_d
        UpdateRK2<<<(point_num / NTH) + 1, NTH>>>(Mhs_d,
                                                  Mhk_d,
                                                  Whs_d,
                                                  Whk_d,
                                                  Wk_d,
                                                  Rhos_d,
                                                  Rhok_d,
                                                  pressures_d,
                                                  pressurek_d,
                                                  func_r_d,
                                                  Altitude_d,
                                                  Altitudeh_d,
                                                  point_num,
                                                  nv);

        BENCH_POINT_I_S_PHY(
            current_step,
            rk,
            "RK2",
            (),
            ("Rhos_d", "Rhok_d", "Mhs_d", "Mhk_d", "Whs_d", "Whk_d", "pressures_d", "pressurek_d"))
    }
    //  Update diagnostic variables.
    cudaDeviceSynchronize();

    BENCH_POINT_I_PHY(
        current_step, "END", (), ("Rho_d", "pressure_d", "Mh_d", "Wh_d", "temperature_d", "W_d"))

    cudaMemcpy(Mh_d, Mhk_d, point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Wh_d, Whk_d, point_num * nvi * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(W_d, Wk_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Rho_d, Rhok_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(pressure_d, pressurek_d, point_num * nv * sizeof(double), cudaMemcpyDeviceToDevice);

    if (phy_modules_execute) phy_modules_dyn_core_loop_end(*this);
}
//END OF THOR!
