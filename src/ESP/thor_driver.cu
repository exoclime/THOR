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

#include <fstream>
#include <iostream>

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
#include "phy/profx_sponge.h"
#include "phy/ultrahot_thermo.h"
#include "reduction_add.h"

#include "diagnostics.h"
#include "phy_modules.h"

__host__ void ESP::Thor(const SimulationSetup& sim, kernel_diagnostics& diag) {
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

    dim3 NBALL((point_num / NTH) + 1, nv, 6); //Number of blocks to execute on all grid points
    dim3 NBALL0((point_num / NTH) + 1, 1, 6); //Number of blocks to execute on all grid points
    dim3 NBRT((point_num / NTH) + 1, 1, 1);
    dim3 NBALL1((point_num / NTH) + 1, nv, 1); //Number of blocks to execute on all grid points

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

    bool energy_equation = false;
    if (thermo_equation == ENERGY) {
        energy_equation = true;
    }

    if (phy_modules_execute)
        phy_modules_dyn_core_loop_init(*this);

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
        if (rk == 0)
            ns_it = 1;
        if (rk == 1)
            ns_it = ns_totali / 2;
        if (rk == 2)
            ns_it = ns_totali;

        if (rk == 0)
            times = timestep / 3.0;
        if (rk == 1)
            times = timestep / ns_totald;
        if (rk == 2)
            times = timestep / ns_totald;

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
                                                            func_r_d,
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
        bool calcT = true;
        if (ultrahot_thermo == VARY_R_CP) {
            check_h = false;
            cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
            update_temperature_Rd_Cp<<<NBALL1, NTH>>>(temperature_d,
                                                      Rd_d,
                                                      Cp_d,
                                                      pressurek_d,
                                                      Rhok_d,
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

        Compute_Temperature_H_Pt_Geff<<<(point_num / NTH) + 1, NTH>>>(temperature_d,
                                                                      pressurek_d,
                                                                      Rhok_d,
                                                                      Mhk_d,
                                                                      h_d,
                                                                      hh_d,
                                                                      pt_d,
                                                                      pth_d,
                                                                      gtil_d,
                                                                      gtilh_d,
                                                                      Whk_d,
                                                                      Wk_d,
                                                                      epotential_d,
                                                                      epotentialh_d,
                                                                      ekinetic_d,
                                                                      ekinetich_d,
                                                                      Etotal_d,
                                                                      sim.P_Ref,
                                                                      sim.Gravit,
                                                                      Cp_d,
                                                                      Rd_d,
                                                                      Altitude_d,
                                                                      Altitudeh_d,
                                                                      point_num,
                                                                      nv,
                                                                      calcT,
                                                                      energy_equation);


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

        if (sim.HyDiff) {
            cudaMemset(diff_d, 0, sizeof(double) * 6 * point_num * nv);

            cudaDeviceSynchronize();
            bool firststep;
            for (int ihyp = 0; ihyp < sim.HyDiffOrder / 2 - 1; ihyp++) {
                if (ihyp == 0)
                    firststep = 1;
                else
                    firststep = 0;
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
                                                  areas_d,
                                                  nvecoa_d,
                                                  nvecti_d,
                                                  nvecte_d,
                                                  func_r_d,
                                                  Kdh4_d,
                                                  Altitude_d,
                                                  sim.A,
                                                  Rd_d,
                                                  Cp_d,
                                                  maps_d,
                                                  nl_region,
                                                  firststep,
                                                  0,
                                                  sim.DeepModel,
                                                  sim.DiffSponge,
                                                  order_diff_sponge,
                                                  Kdh2_d,
                                                  boundary_flux_d,
                                                  energy_equation,
                                                  sim.HyDiffOrder);

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
                                                   areas_d,
                                                   nvecoa_d,
                                                   nvecti_d,
                                                   nvecte_d,
                                                   Kdh4_d,
                                                   Altitude_d,
                                                   Altitudeh_d,
                                                   sim.A,
                                                   Rd_d,
                                                   Cp_d,
                                                   point_local_d,
                                                   point_num,
                                                   firststep,
                                                   0,
                                                   sim.DeepModel,
                                                   sim.DiffSponge,
                                                   order_diff_sponge,
                                                   Kdh2_d,
                                                   boundary_flux_d,
                                                   energy_equation,
                                                   sim.HyDiffOrder);
                cudaDeviceSynchronize();
            }

            BENCH_POINT_I_S(current_step, rk, "Diffusion_Op1", (), ("diff_d"))
            //  ("diffmh_d", "diffw_d", "diffrh_d", "diffpr_d", "diff_d"))


            cudaMemset(diffrh_d, 0, sizeof(double) * point_num * nv);
            cudaMemset(boundary_flux_d, 0, sizeof(double) * 6 * point_num * nv);
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
                                              areas_d,
                                              nvecoa_d,
                                              nvecti_d,
                                              nvecte_d,
                                              func_r_d,
                                              Kdh4_d,
                                              Altitude_d,
                                              sim.A,
                                              Rd_d,
                                              Cp_d,
                                              maps_d,
                                              nl_region,
                                              0,
                                              1,
                                              sim.DeepModel,
                                              sim.DiffSponge,
                                              order_diff_sponge,
                                              Kdh2_d,
                                              boundary_flux_d,
                                              energy_equation,
                                              sim.HyDiffOrder);
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
                                               areas_d,
                                               nvecoa_d,
                                               nvecti_d,
                                               nvecte_d,
                                               Kdh4_d,
                                               Altitude_d,
                                               Altitudeh_d,
                                               sim.A,
                                               Rd_d,
                                               Cp_d,
                                               point_local_d,
                                               point_num,
                                               0,
                                               1,
                                               sim.DeepModel,
                                               sim.DiffSponge,
                                               order_diff_sponge,
                                               Kdh2_d,
                                               boundary_flux_d,
                                               energy_equation,
                                               sim.HyDiffOrder);

            cudaDeviceSynchronize();

            BENCH_POINT_I_S(current_step,
                            rk,
                            "Diffusion_Op2",
                            (),
                            ("diffmh_d", "diffw_d", "diffrh_d", "diffpr_d", "diff_d"))

            if (sim.VertHyDiff) { //sixth-order vertical hyperdiffusion
                cudaMemset(diff_d, 0, sizeof(double) * 6 * point_num * nv);
                cudaMemset(diff2_d, 0, sizeof(double) * 6 * point_num * nv);

                vertical_diff<<<NBALL0, NTH>>>(diffmv_d,  //
                                               diffwv_d,  //
                                               diffrv_d,  //
                                               diffprv_d, //
                                               diff_d,
                                               diff2_d,
                                               Mhk_d,         //
                                               Rhok_d,        //
                                               temperature_d, //
                                               Wk_d,          //
                                               Kdv6_d,        //
                                               Altitude_d,    //
                                               Altitudeh_d,
                                               sim.A,  //
                                               sim.Rd, //
                                               sim.VertHyDiffOrder,
                                               sim.DeepModel,
                                               point_num,
                                               nv);

                // vertical_diff_joao<40><<<NBALL0, NTH>>>(diffmv_d,      //
                //                                         diffwv_d,      //
                //                                         diffrv_d,      //
                //                                         diffprv_d,     //
                //                                         Mhk_d,         //
                //                                         Rhok_d,        //
                //                                         temperature_d, //
                //                                         Wk_d,          //
                //                                         Kdv6_d,        //
                //                                         Altitude_d,    //
                //                                         Altitudeh_d,
                //                                         sim.A,  //
                //                                         sim.Rd, //
                //                                         sim.DeepModel,
                //                                         point_num);
                cudaDeviceSynchronize();
            }

            Correct_Horizontal<<<NBALL1, NTH>>>(diffmh_d, diffmv_d, func_r_d, point_num);

            cudaDeviceSynchronize();

            BENCH_POINT_I_S(
                current_step, rk, "Correct_Horizontal", (), ("diffmh_d", "diffmv_d", "func_r_d"))
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
        
        bool increased_damping_for_100_days = 0;
        if (increased_damping_for_100_days == 1) {
            

            if (current_step==1) {

                // pascicode
    
                
                double Diffc_intitial_factor = 0.1;
                double DivDampc_intitial_factor = 1.0;
                double Diffc_v_intitial_factor = 1.0;

                printf("increased Diffc in the first 100 time steps by a factor of %e for PARMENTIER \n",Diffc_intitial_factor);

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
                    Kdh4_h[lev] = (Diffc_intitial_factor*sim.Diffc) * pow(dbar, 1.0 * sim.HyDiffOrder)
                                / timestep; // * Altitude_h[lev]/sim.Top_altitude;
                    Kdhz_h[lev] =
                        (DivDampc_intitial_factor*sim.DivDampc) * pow(dbar, 4.) / timestep; // * Altitude_h[lev]/sim.Top_altitude;
                    if (sim.DiffSponge) {
                        double n = Altitude_h[lev] / sim.Top_altitude;
                        if (n > ns_diff_sponge) {
                            ksponge = Dv_sponge
                                    * pow(sin(0.5 * M_PI * (n - ns_diff_sponge) / (1.0 - ns_diff_sponge)), 2);
                        }
                        else {
                            ksponge = 0;
                        }
                        if (order_diff_sponge == 2) {
                            Kdh2_h[lev] = ksponge * pow(dbar, 2.) / timestep;
                        }
                        else if (order_diff_sponge == 4) {
                            Kdh4_h[lev] += ksponge * pow(dbar, 4.) / timestep;
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
                    double dz   = Altitudeh_h[lev + 1] - Altitudeh_h[lev];
                    Kdv6_h[lev] = Diffc_v_intitial_factor * sim.Diffc_v * pow(dz, 1.0 * sim.VertHyDiffOrder) / timestep;
                    Kdvz_h[lev] = 0.0; //not used (yet? perhaps in future)
                }
    
                cudaMemcpy(Kdhz_d, Kdhz_h, nv * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(Kdh4_d, Kdh4_h, nv * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(Kdv6_d, Kdv6_h, nv * sizeof(double), cudaMemcpyHostToDevice);
    
            }
            if (current_step==100) {
                
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
                    Kdh4_h[lev] = (sim.Diffc) * pow(dbar, 1.0 * sim.HyDiffOrder)
                                / timestep; // * Altitude_h[lev]/sim.Top_altitude;
                    Kdhz_h[lev] =
                        (sim.DivDampc) * pow(dbar, 4.) / timestep; // * Altitude_h[lev]/sim.Top_altitude;
                    if (sim.DiffSponge) {
                        double n = Altitude_h[lev] / sim.Top_altitude;
                        if (n > ns_diff_sponge) {
                            ksponge = Dv_sponge
                                    * pow(sin(0.5 * M_PI * (n - ns_diff_sponge) / (1.0 - ns_diff_sponge)), 2);
                        }
                        else {
                            ksponge = 0;
                        }
                        if (order_diff_sponge == 2) {
                            Kdh2_h[lev] = ksponge * pow(dbar, 2.) / timestep;
                        }
                        else if (order_diff_sponge == 4) {
                            Kdh4_h[lev] += ksponge * pow(dbar, 4.) / timestep;
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
                    double dz   = Altitudeh_h[lev + 1] - Altitudeh_h[lev];
                    Kdv6_h[lev] = sim.Diffc_v * pow(dz, 1.0 * sim.VertHyDiffOrder) / timestep;
                    Kdvz_h[lev] = 0.0; //not used (yet? perhaps in future)
                }
    
                cudaMemcpy(Kdhz_d, Kdhz_h, nv * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(Kdh4_d, Kdh4_h, nv * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(Kdv6_d, Kdv6_h, nv * sizeof(double), cudaMemcpyHostToDevice);
            }

        }
        
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
                                        areas_d,
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
                                         areas_d,
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
                                        areas_d,
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
                                         areas_d,
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

        if (sim.RayleighSponge == true && raysp_calc_mode != IMP) {
            cudaMemset(profx_dMh_d, 0, sizeof(double) * 3 * point_num * nv);
            cudaMemset(profx_dWh_d, 0, sizeof(double) * point_num * nvi);
            cudaMemset(profx_dW_d, 0, sizeof(double) * point_num * nv);
            if (rk == 0 || raysp_calc_mode == EXP3) {
                // dim3 NBT((point_num / NTH) + 1, nv, 1);

                if (damp_uv_to_mean) {
                    cudaMemset(utmp, 0, sizeof(double) * nlat_bins * nv * max_count);
                    cudaMemset(vtmp, 0, sizeof(double) * nlat_bins * nv * max_count);

                    zonal_uv<<<NB, NTH>>>(Mhk_d,
                                          Rhok_d,
                                          zonal_mean_tab_d,
                                          lonlat_d,
                                          point_num,
                                          utmp,
                                          vtmp,
                                          max_count);

                    cudaDeviceSynchronize();
                }
                if (damp_w_to_mean) {
                    cudaMemset(wtmp, 0, sizeof(double) * nlat_bins * nv * max_count);
                    zonal_w<<<NB, NTH>>>(
                        Wk_d, Rhok_d, zonal_mean_tab_d, point_num, wtmp, max_count);

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
                    cudaMemcpy(vbar_d,
                               vbar_h,
                               3 * nlat_bins * nv * sizeof(double),
                               cudaMemcpyHostToDevice);
                }

                if (sim.RayleighSpongeT) {
                    zonal_temp<<<NB, NTH>>>(pressurek_d,
                                            Rhok_d,
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
                    cudaMemcpy(
                        Tbar_d, Tbar_h, nlat_bins * nv * sizeof(double), cudaMemcpyHostToDevice);
                }
            }
            double Rv_fac = 1;
            if (shrink_sponge == true) {
                if (current_step * timestep >= t_shrink * timestep) {
                    double shrink_scale = timestep * 1000;
                    Rv_fac = exp(-(current_step * timestep - t_shrink * timestep) / shrink_scale);
                }
            }


            sponge_layer<<<NBRT, NTH>>>(Mhk_d,
                                        Rhok_d,
                                        Wk_d,
                                        Whk_d,
                                        pressurek_d,
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
                                        false,
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

            BENCH_POINT_I_S(current_step,
                            rk,
                            "phy_Sponge",
                            (),
                            ("Rhok_d", "pressurek_d", "Mhk_d", "Whk_d", "temperature_d", "Wk_d"))

            cudaDeviceSynchronize();
        }

        BENCH_POINT_I_S(current_step,
                        rk,
                        "phy_before_slow_modews",
                        (),
                        ("Mhk_d",
                         "Whk_d",
                         "Rhok_d",
                         "Adv_d",
                         "DivM_d",
                         "diffmh_d",
                         "diffw_d",
                         "diffrh_d",
                         "diffpr_d",
                         "diffmv_d",
                         "diffwv_d",
                         "diffrv_d",
                         "diffprv_d",
                         "pressurek_d",
                         "h_d",
                         "hh_d",
                         "gtil_d",
                         "Altitude_d",
                         "Altitudeh_d",
                         "Qheat_d"))

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
                                               Cp_d,
                                               Rd_d,
                                               func_r_d,
                                               maps_d,
                                               nl_region,
                                               sim.DeepModel,
                                               sim.NonHydro,
                                               profx_dMh_d,
                                               profx_dWh_d,
                                               profx_Qheat_d,
                                               energy_equation);
        cudaDeviceSynchronize();
        BENCH_POINT_I_S(current_step,
                        rk,
                        "Compute_Slow_Modes",
                        (),
                        ("SlowMh_d", "SlowWh_d", "SlowRho_d", "Slowpressure_d"))
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
                                              Cp_d,
                                              Rd_d,
                                              func_r_d,
                                              point_local_d,
                                              nv,
                                              point_num,
                                              sim.DeepModel,
                                              sim.NonHydro,
                                              profx_dMh_d,
                                              profx_dWh_d,
                                              profx_Qheat_d,
                                              energy_equation);
        cudaDeviceSynchronize();

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
        //      SMALL-STEP
        for (int ns = 0; ns < ns_it; ns++) {
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
                                            areas_d,
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
                                             areas_d,
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
                                            areas_d,
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
                                             areas_d,
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
                                                          Cp_d,
                                                          Rd_d,
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
                                                         Cp_d,
                                                         Rd_d,
                                                         sim.A,
                                                         point_local_d,
                                                         point_num,
                                                         nv,
                                                         sim.DeepModel);

            BENCH_POINT_I_SS(
                current_step, rk, ns, "Prepare_Implicit_Vertical", (), ("Sp_d", "Sd_d"))

            cudaDeviceSynchronize();
#if defined(DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_DIAG_DOM) \
    || defined(DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_RESULT)
            // reset diagnostics flag
            // diagnostics array is reset inside kernel
            // could also be reset here with diag.reset();
            diag.reset_flag();
#endif // DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_DIAG_DOM

            // Updates: Whs_d, Ws_d
            Vertical_Eq<<<(point_num / num_th_vertical_eq) + 1,
                          num_th_vertical_eq,
                          2 * num_th_vertical_eq * nvi * sizeof(double)>>>(
                Whs_d,
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
                Cp_d,
                Rd_d,
                times,
                sim.Gravit,
                Altitude_d,
                Altitudeh_d,
                sim.A,
                sim.NonHydro,
                point_num,
                nv,
                nvi,
                sim.DeepModel,
                *diag.diagnostics_global_flag, // Pass the diagnostics arrays device pointers
                *diag.diagnostics);

            //          Pressure and density equations.
            cudaDeviceSynchronize();
            BENCH_POINT_I_SS_PHY(current_step,
                                 rk,
                                 ns,
                                 "Vertical_Eq",
                                 (),
                                 ("Whs_d", "Ws_d", "pressures_d", "h_d", "hh_d", "Rhos_d"));

#if defined(DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_DIAG_DOM) \
    || defined(DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_RESULT)
            if (diag.check_flag()) {
                // dump the data
                diag.dump_data("vertical_eq", current_step, rk, ns, *this, point_num, nvi);
                // warn, but do not quit.
                unsigned int flag = diag.get_flag();
                // Check for global flags
                if ((flag & THOMAS_NOT_DD) == THOMAS_NOT_DD)
                    log::printf(
                        "Non diagonaly dominant matrix for thomas algorithm in Vertical_Eq\n");
                if ((flag & THOMAS_BAD_SOLUTION) == THOMAS_BAD_SOLUTION)
                    log::printf("Bad thomas algorithm solution in Vertical_Eq\n");
            }
#endif // DIAG_CHECK_THOR_VERTICAL_INT_THOMAS_DIAG_DOM

            cudaError_t err = cudaGetLastError();

            // Check device query
            if (err != cudaSuccess) {
                log::printf("[%s:%d] CUDA error check reports error: %s\n",
                            __FILE__,
                            __LINE__,
                            cudaGetErrorString(err));
                exit(EXIT_FAILURE);
            }

            // update the physics modules in fast mode
            if (phy_modules_execute)
                phy_modules_dyn_core_loop_fast_modes(*this, sim, current_step, times);


            BENCH_POINT_I_SS_PHY(current_step,
                                 rk,
                                 ns,
                                 "Phy_mod_fast_mode",
                                 (),
                                 ("Whs_d", "Ws_d", "pressures_d", "h_d", "hh_d", "Rhos_d"))

#if defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX) || defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN)
            diag.reset_flag();
#endif
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
                                                     epotential_d,
                                                     epotentialh_d,
                                                     ekinetic_d,
                                                     ekinetich_d,
                                                     Etotal_d,
                                                     h_d,
                                                     hh_d,
                                                     SlowRho_d,
                                                     profx_Qheat_d,
                                                     diffpr_d,
                                                     diffprv_d,
                                                     div_d,
                                                     Altitude_d,
                                                     Altitudeh_d,
                                                     Cp_d,
                                                     Rd_d,
                                                     sim.A,
                                                     sim.P_Ref,
                                                     sim.Gravit,
                                                     times,
                                                     maps_d,
                                                     nl_region,
                                                     sim.DeepModel,
                                                     energy_equation,
                                                     *diag.diagnostics_global_flag,
                                                     *diag.diagnostics);


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
                                                    epotential_d,
                                                    epotentialh_d,
                                                    ekinetic_d,
                                                    ekinetich_d,
                                                    Etotal_d,
                                                    h_d,
                                                    hh_d,
                                                    SlowRho_d,
                                                    profx_Qheat_d,
                                                    diffpr_d,
                                                    diffprv_d,
                                                    div_d,
                                                    Altitude_d,
                                                    Altitudeh_d,
                                                    Cp_d,
                                                    Rd_d,
                                                    sim.A,
                                                    sim.P_Ref,
                                                    sim.Gravit,
                                                    times,
                                                    point_local_d,
                                                    point_num,
                                                    nv,
                                                    sim.DeepModel,
                                                    energy_equation,
                                                    *diag.diagnostics_global_flag,
                                                    *diag.diagnostics);

            BENCH_POINT_I_SS(
                current_step, rk, ns, "Density_Pressure_Eqs_Poles", (), ("pressures_d", "Rhos_d"))

#if defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_AUX) || defined(DIAG_CHECK_DENSITY_PRESSURE_EQ_P_NAN)
            if (diag.check_flag()) {
                diag.dump_data("density_pressure_eq", current_step, rk, ns, *this, point_num, nv);
                log::printf("NaN value or Negative AUX in Density_Pressure_Eqs, EXITING\n");
                exit(EXIT_FAILURE);
            }
#endif
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
        cudaDeviceSynchronize();

        BENCH_POINT_I_S_PHY(
            current_step,
            rk,
            "RK2",
            (),
            ("Rhos_d", "Rhok_d", "Mhs_d", "Mhk_d", "Whs_d", "Whk_d", "pressures_d", "pressurek_d"))

        if (ultrahot_thermo == VARY_R_CP) {
            check_h = false;
            cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
            update_temperature_Rd_Cp<<<NBALL1, NTH>>>(temperature_d,
                                                      Rd_d,
                                                      Cp_d,
                                                      pressurek_d,
                                                      Rhok_d,
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
        }
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

    if (phy_modules_execute)
        phy_modules_dyn_core_loop_end(*this);
}
//END OF THOR!
