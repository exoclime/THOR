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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "../headers/esp.h"                   // Global parameters.

#include "../headers/dyn/thor_adv_cor.h"      // Advection term.
#include "../headers/dyn/thor_auxiliary.h"    // Temperature, interal energy, potential tempareture and effective gravity.
#include "../headers/dyn/thor_diff.h"         // Hyper-diffusion.
#include "../headers/dyn/thor_div.h"          // Divergence damping.
#include "../headers/dyn/thor_fastmodes.h"    // Fast terms.
#include "../headers/dyn/thor_slowmodes.h"    // Slow terms.
#include "../headers/dyn/thor_vertical_int.h" // Vertical momentum.


#include "binary_test.h"

__host__ void ESP::Thor(double timestep_dyn, // Large timestep.
                        bool   HyDiff      , // Turn on/off hyper-diffusion.
                        bool   DivDampP    , // Turn on/off divergence damping.
                        double Omega       , // Rotation rate.
                        double Cp          , // Heat capaciry.
                        double Rd          , // Gas constant (atmosphere).
                        double Mmol        , // Molecular mass (atmosphere).
                        double mu          , // Mass unit.
                        double kb          , // Boltzmann constant.
                        double P_Ref       , // Averaged pressure surface.
                        double Gravit      , // Gravity.
                        double A           , // Planet radius.
                        bool NonHydro      , // Turn on/off non-hydrostatic.
                        bool DeepModel     ){// Turn on/off deep atmosphere.
//
//  Number of threads per block.
    const int NTH = 256;

    // Vertical Eq only works on vertical stack of data, can run independently, only uses shared
    // memory for intermediate data that is not shared with neighbours.
    // Need to set the block size so that the internal arrays fit in shared memory for each block
    const int num_th_vertical_eq = 32;
//  Specify the block sizes.
    const int LN = 16;               // Size of the inner region side.
    dim3 NT(nl_region, nl_region, 1);// Number of threads in a block.
    dim3 NB(nr, nv , 1);             // Number of blocks.
    dim3 NBD(nr, nv, 6);             // Number of blocks in the diffusion routine.
    dim3 NBDP(2, nv, 6);             // Number of blocks in the diffusion routine. (POLES)
    dim3 NBP(2, nv, 1);              // Number of blocks. (POLES)

//  Number of Small steps
    double ns_totald = 6;            // Maximum number of small steps in a large step (double ).
    int    ns_totali = 6;            // Maximum number of small steps in a large step (integer).
    int    ns_it        ;            // Number of small steps in each large step.
    double times        ;            // Sub-timestep.

//  Initialize local variables used for the time integration.
    cudaDeviceSynchronize();
    cudaMemcpy(Mhk_d       , Mh_d       , point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Whk_d       , Wh_d       , point_num * nvi*     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Wk_d        , W_d        , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Rhok_d      , Rho_d      , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(pressurek_d , pressure_d , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);


    cudaMemset(Mhs_d       , 0, sizeof(double) * 3*point_num * nv);
    cudaMemset(Rhos_d      , 0, sizeof(double) * point_num * nv)  ;
    cudaMemset(Whs_d       , 0, sizeof(double) * point_num * nvi) ;
    cudaMemset(Ws_d        , 0, sizeof(double) * point_num * nv)  ;
    cudaMemset(pressures_d , 0, sizeof(double) * point_num * nv)  ;

    USE_BENCHMARK()

    BENCH_POINT_I(*this, current_step, "thor_init")

//  Loop for large time integration.
    for(int rk = 0; rk < 3; rk++){
//      Local variables to define the length (times) and the number of the small steps (ns_it).
        if(rk == 0) ns_it = 1           ;
        if(rk == 1) ns_it = ns_totali/2 ;
        if(rk == 2) ns_it = ns_totali   ;

        if(rk == 0) times = timestep_dyn/3.0      ;
        if(rk == 1) times = timestep_dyn/ns_totald;
        if(rk == 2) times = timestep_dyn/ns_totald;

        // initialise some memory

//
//      Compute advection and coriolis terms.
        cudaMemset(Adv_d, 0, sizeof(double) * 3 * point_num * nv); // Sets every value of Adv_d to
                                                                   // zero.
        cudaDeviceSynchronize();


        // Updates: Adv_d, v_d
        Compute_Advec_Cori1<LN,LN><<<NB,NT>>>(Adv_d      ,
                                              v_d        ,
                                              Mhk_d      ,
                                              div_d      ,
                                              Wk_d       ,
                                              Rhok_d     ,
                                              Altitude_d ,
                                              A          ,
                                              func_r_d   ,
                                              maps_d     ,
                                              nl_region  ,
                                              DeepModel);
        // Updates: Adv_d, v_d
        Compute_Advec_Cori_Poles<6><<<2, 1>>>(Adv_d        ,
                                              v_d          ,
                                              Mhk_d        ,
                                              div_d        ,
                                              Wk_d         ,
                                              Rhok_d       ,
                                              Altitude_d   ,
                                              A            ,
                                              func_r_d     ,
                                              point_local_d,
                                              point_num    ,
                                              nv           ,
                                              DeepModel   );

        cudaDeviceSynchronize();
        // Updates: Adv_d
        Compute_Advec_Cori2 <<<(point_num / NTH) + 1, NTH >>>(Adv_d      ,
                                                              v_d        ,
                                                              Whk_d      ,
                                                              Rhok_d     ,
                                                              Altitude_d ,
                                                              Altitudeh_d,
                                                              Omega      ,
                                                              A          ,
                                                              nv         ,
                                                              point_num  ,
                                                              DeepModel  );

        check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_thor<<< 16, NTH >>>(Adv_d, nv, 3*point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(v_d, nv, 3*point_num, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
        if(check_h){
           printf("\n\n Error in NAN check after Thor:Compute_Advec_Cori!\n");
           exit(EXIT_FAILURE);
        }
//
//      Computes temperature, internal energy, potential temperature and effective gravity.
        cudaDeviceSynchronize();

        BENCH_POINT_I_S(*this, current_step, rk, "Compute_Advec_Cori")

        // Updates: temperature_d, h_d, hh_d, pt_d, pth_d, gtil_d, gtilh_d
        Compute_Temperature_H_Pt_Geff <<< (point_num / NTH) + 1, NTH >>> (temperature_d,
                                                                          pressurek_d  ,
                                                                          Rhok_d       ,
                                                                          h_d          ,
                                                                          hh_d         ,
                                                                          pt_d         ,
                                                                          pth_d        ,
                                                                          gtil_d       ,
                                                                          gtilh_d      ,
                                                                          Whk_d        ,
                                                                          P_Ref        ,
                                                                          Gravit       ,
                                                                          Cp           ,
                                                                          Rd           ,
                                                                          Altitude_d   ,
                                                                          Altitudeh_d  ,
                                                                          point_num    ,
                                                                          nv           );

        check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_thor<<< 16, NTH >>>(temperature_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(h_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(hh_d, nv+1, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(pt_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(pth_d, nv+1, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(gtil_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(gtilh_d, nv+1, point_num, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
        if(check_h){
           printf("\n\n Error in NAN check after Thor:Compute_Temp!\n");
           exit(EXIT_FAILURE);
        }
//      Initializes slow terms.
        cudaDeviceSynchronize();
        cudaMemset(SlowMh_d      , 0, sizeof(double) * 3 * point_num * nv);
        cudaMemset(SlowWh_d      , 0, sizeof(double) *     point_num * nvi);
        cudaMemset(SlowRho_d     , 0, sizeof(double) *     point_num * nv);
        cudaMemset(Slowpressure_d, 0, sizeof(double) *     point_num * nv);
//
//      Hyper-Diffusion.
        if (HyDiff){
            cudaMemset(diff_d, 0, sizeof(double) *   6 * point_num * nv);
            cudaDeviceSynchronize();
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op <LN,LN> <<< NBD,NT >>>(diffmh_d     ,
                                                diffw_d      ,
                                                diffrh_d     ,
                                                diffpr_d     ,
                                                diff_d       ,
                                                Mhk_d        ,
                                                Rhok_d       ,
                                                temperature_d,
                                                Wk_d         ,
                                                areasTr_d    ,
                                                nvecoa_d     ,
                                                nvecti_d     ,
                                                nvecte_d     ,
                                                func_r_d     ,
                                                Kdh4_d       ,
                                                Altitude_d   ,
                                                A            ,
                                                Rd           ,
                                                maps_d       ,
                                                nl_region    ,
                                                0            ,
                                                DeepModel    );
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op_Poles <5><<<NBDP,1 >>>(diffmh_d     ,
                                                 diffw_d      ,
                                                diffrh_d     ,
                                                diffpr_d     ,
                                                diff_d       ,
                                                Mhk_d        ,
                                                Rhok_d       ,
                                                temperature_d,
                                                Wk_d         ,
                                                func_r_d     ,
                                                areasTr_d    ,
                                                nvecoa_d     ,
                                                nvecti_d     ,
                                                nvecte_d     ,
                                                Kdh4_d       ,
                                                Altitude_d   ,
                                                Altitudeh_d  ,
                                                A            ,
                                                Rd           ,
                                                point_local_d,
                                                point_num    ,
                                                0            ,
                                                DeepModel    );
            cudaDeviceSynchronize();
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op <LN,LN><<<  NBD,NT >>>(diffmh_d     ,
                                                diffw_d      ,
                                                diffrh_d     ,
                                                  diffpr_d     ,
                                                diff_d       ,
                                                Mhk_d        ,
                                                Rhok_d       ,
                                                temperature_d,
                                                Wk_d         ,
                                                areasTr_d    ,
                                                nvecoa_d     ,
                                                nvecti_d     ,
                                                nvecte_d     ,
                                                func_r_d     ,
                                                Kdh4_d       ,
                                                Altitude_d   ,
                                                A            ,
                                                Rd           ,
                                                maps_d       ,
                                                nl_region    ,
                                                1            ,
                                                DeepModel    );
            //Updates: diffmh_d, diffw_d, diffrh_d, diffpr_d, diff_d
            Diffusion_Op_Poles <5><<<NBDP, 1>>>(diffmh_d     ,
                                                 diffw_d      ,
                                                diffrh_d     ,
                                                diffpr_d     ,
                                                diff_d       ,
                                                Mhk_d        ,
                                                Rhok_d       ,
                                                temperature_d,
                                                Wk_d         ,
                                                func_r_d     ,
                                                areasTr_d    ,
                                                nvecoa_d     ,
                                                nvecti_d     ,
                                                nvecte_d     ,
                                                Kdh4_d       ,
                                                Altitude_d   ,
                                                Altitudeh_d  ,
                                                A            ,
                                                Rd           ,
                                                point_local_d,
                                                point_num    ,
                                                1            ,
                                                DeepModel    );
            check_h = false;
            cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
            isnan_check_thor<<< 16, NTH >>>(diffmh_d, nv, 3*point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(diffw_d, nv, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(diffrh_d, nv, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(diffpr_d, nv, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(diff_d, nv, 6*point_num, check_d);
            cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
            if(check_h){
               printf("\n\n Error in NAN check after Thor:Diffusion_Op!\n");
               exit(EXIT_FAILURE);
            }
        }


//
//      Divergence damping
        cudaMemset(DivM_d   , 0, sizeof(double)*point_num * 3 * nv);
        cudaMemset(divg_Mh_d, 0, sizeof(double)*point_num * 3 * nv);
        if (DivDampP) {
            cudaDeviceSynchronize();
            // Updates: DivM_d, divg_Mh_d
            DivM_Op <LN,LN> <<< NB, NT >>>(DivM_d     ,
                                           divg_Mh_d  ,
                                           Mhk_d      ,
                                           Whk_d      ,
                                           Kdhz_d     ,
                                           areasTr_d  ,
                                           nvecoa_d   ,
                                           nvecti_d   ,
                                           nvecte_d   ,
                                           func_r_d   ,
                                           Altitudeh_d,
                                           Altitude_d ,
                                           A          ,
                                           maps_d     ,
                                           nl_region  ,
                                           0          ,
                                           DeepModel);
            // Updates: DivM_d, divg_Mh_d
            DivM_Op_Poles<5><<<NBP,1>>>(DivM_d       ,
                                        divg_Mh_d    ,
                                        Mhk_d        ,
                                        Whk_d        ,
                                        Kdhz_d       ,
                                        areasTr_d    ,
                                        nvecoa_d     ,
                                        nvecti_d     ,
                                        nvecte_d     ,
                                        func_r_d     ,
                                        Altitudeh_d  ,
                                        Altitude_d   ,
                                        A            ,
                                        point_local_d,
                                        point_num    ,
                                        0            ,
                                        DeepModel    );

            cudaDeviceSynchronize();
            // Updates: DivM_d, divg_Mh_d
            DivM_Op <LN,LN> <<< NB, NT >>>(DivM_d     ,
                                           divg_Mh_d  ,
                                           Mhk_d      ,
                                           Whk_d      ,
                                           Kdhz_d     ,
                                           areasTr_d  ,
                                           nvecoa_d   ,
                                           nvecti_d   ,
                                           nvecte_d   ,
                                           func_r_d   ,
                                           Altitudeh_d,
                                           Altitude_d ,
                                           A          ,
                                           maps_d     ,
                                           nl_region  ,
                                           1          ,
                                           DeepModel);
            // Updates: DivM_d, divg_Mh_d
            DivM_Op_Poles<5><<<NBP,1>>>(DivM_d       ,
                                        divg_Mh_d    ,
                                        Mhk_d        ,
                                        Whk_d        ,
                                        Kdhz_d       ,
                                        areasTr_d    ,
                                        nvecoa_d     ,
                                        nvecti_d     ,
                                        nvecte_d     ,
                                        func_r_d     ,
                                        Altitudeh_d  ,
                                        Altitude_d   ,
                                        A            ,
                                        point_local_d,
                                        point_num    ,
                                        1            ,
                                        DeepModel    );
            check_h = false;
            cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
            isnan_check_thor<<< 16, NTH >>>(DivM_d, nv, 3*point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(divg_Mh_d, nv, 3*point_num, check_d);
            cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
            if(check_h){
               printf("\n\n Error in NAN check after Thor:DivM_Op slow!\n");
               exit(EXIT_FAILURE);
            }
        }

//
//      Slow Modes
        cudaDeviceSynchronize();
        // Updates: SlowMh_d, SlowWh_d, SlowRho_d, Slowpressure_d
        Compute_Slow_Modes <LN,LN>  <<<NB, NT >>> (SlowMh_d      ,
                                                   SlowWh_d      ,
                                                   SlowRho_d     ,
                                                   Slowpressure_d,
                                                   Mhk_d         ,
                                                   Whk_d         ,
                                                   Rhok_d        ,
                                                   Adv_d         ,
                                                   DivM_d        ,
                                                   diffmh_d      ,
                                                   diffw_d       ,
                                                   diffrh_d      ,
                                                   diffpr_d      ,
                                                   pressurek_d   ,
                                                   h_d           ,
                                                   hh_d          ,
                                                   gtil_d        ,
                                                   grad_d        ,
                                                   div_d         ,
                                                   Altitude_d    ,
                                                   Altitudeh_d   ,
                                                   A             ,
                                                   Gravit        ,
                                                   Cp            ,
                                                   Rd            ,
                                                   func_r_d      ,
                                                   maps_d        ,
                                                   nl_region     ,
                                                   DeepModel     ,
                                                   NonHydro      );
        cudaDeviceSynchronize();
        // Updates: SlowMh_d, SlowWh_d, SlowRho_d, Slowpressure_d
        Compute_Slow_Modes_Poles<6><<<2,1>>> (SlowMh_d,
                                              SlowWh_d      ,
                                              SlowRho_d     ,
                                              Slowpressure_d,
                                              Mhk_d         ,
                                              Whk_d         ,
                                              Rhok_d        ,
                                              Adv_d         ,
                                              DivM_d        ,
                                              diffmh_d      ,
                                              diffw_d       ,
                                              diffrh_d      ,
                                              diffpr_d      ,
                                              pressurek_d   ,
                                              h_d           ,
                                              hh_d          ,
                                              gtil_d        ,
                                              grad_d        ,
                                              div_d         ,
                                              Altitude_d    ,
                                              Altitudeh_d   ,
                                              A             ,
                                              Gravit        ,
                                              Cp            ,
                                              Rd            ,
                                              func_r_d      ,
                                              point_local_d ,
                                              nv            ,
                                              point_num     ,
                                              DeepModel     ,
                                              NonHydro      );

//      Updates or initializes deviations.
        check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_thor<<< 16, NTH >>>(SlowMh_d, nv, 3*point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(SlowWh_d, nvi, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(SlowRho_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(Slowpressure_d, nv, point_num, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
        if(check_h){
           printf("\n\n Error in NAN check after Thor:SlowModes!\n");
           exit(EXIT_FAILURE);
        }

        if (rk > 0)
        {
            cudaDeviceSynchronize();
            // Updates: Mhs_d, Whs_d, Ws_d, Rhos_d, pressures_d
            UpdateRK <<< (point_num / NTH) + 1, NTH >>> (Mhs_d      ,
                                                         Mhk_d      ,
                                                         Mh_d       ,
                                                         Whs_d      ,
                                                         Whk_d      ,
                                                         Wh_d       ,
                                                         Ws_d       ,
                                                         Rhos_d     ,
                                                         Rhok_d     ,
                                                         Rho_d      ,
                                                         pressures_d,
                                                         pressurek_d,
                                                         pressure_d ,
                                                         func_r_d   ,
                                                         Altitude_d ,
                                                         Altitudeh_d,
                                                         point_num  ,
                                                         nv         );
            check_h = false;
            cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
            isnan_check_thor<<< 16, NTH >>>(Mhs_d, nv, 3*point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(Whs_d, nvi, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(Ws_d, nv, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(Rhos_d, nv, point_num, check_d);
            isnan_check_thor<<< 16, NTH >>>(pressures_d, nv, point_num, check_d);
            cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
            if(check_h){
               printf("\n\n Error in NAN check after Thor:UpdateRK!\n");
               exit(EXIT_FAILURE);
            }
        }

//
//      SMALL-STEPS
        for (int ns = 0; ns < ns_it; ns++){
            cudaMemset(DivM_d   , 0, sizeof(double)*point_num * 3 * nv);
            cudaMemset(divg_Mh_d, 0, sizeof(double)*point_num * 3 * nv);
            if (DivDampP) {
                cudaDeviceSynchronize();
                // Updates: DivM_d, divg_Mh_d
                DivM_Op <LN,LN><<< NB, NT >>>(DivM_d     ,
                                                divg_Mh_d  ,
                                               Mhs_d      ,
                                               Whs_d      ,
                                               Kdhz_d     ,
                                               areasTr_d  ,
                                               nvecoa_d   ,
                                               nvecti_d   ,
                                               nvecte_d   ,
                                               func_r_d   ,
                                               Altitudeh_d,
                                               Altitude_d ,
                                               A          ,
                                               maps_d     ,
                                               nl_region  ,
                                               0          ,
                                               DeepModel);
                // Updates: DivM_d, divg_Mh_d
                DivM_Op_Poles<5><<<NBP,1>>>(DivM_d       ,
                                            divg_Mh_d    ,
                                            Mhs_d        ,
                                            Whs_d        ,
                                            Kdhz_d       ,
                                            areasTr_d    ,
                                            nvecoa_d     ,
                                            nvecti_d     ,
                                            nvecte_d     ,
                                            func_r_d     ,
                                            Altitudeh_d  ,
                                            Altitude_d   ,
                                            A            ,
                                            point_local_d,
                                            point_num    ,
                                            0            ,
                                            DeepModel    );

                cudaDeviceSynchronize();
                // Updates: DivM_d, divg_Mh_d
                DivM_Op <LN, LN> <<< NB, NT>>>(DivM_d,
                                                divg_Mh_d  ,
                                               Mhs_d      ,
                                               Whs_d      ,
                                               Kdhz_d     ,
                                               areasTr_d  ,
                                               nvecoa_d   ,
                                               nvecti_d   ,
                                               nvecte_d   ,
                                               func_r_d   ,
                                               Altitudeh_d,
                                               Altitude_d ,
                                               A          ,
                                               maps_d     ,
                                               nl_region  ,
                                               1          ,
                                               DeepModel);
                // Updates: DivM_d, divg_Mh_d
                DivM_Op_Poles<5><<<NBP,1>>>(DivM_d       ,
                                            divg_Mh_d    ,
                                            Mhs_d        ,
                                            Whs_d        ,
                                            Kdhz_d       ,
                                            areasTr_d    ,
                                            nvecoa_d     ,
                                            nvecti_d     ,
                                            nvecte_d     ,
                                            func_r_d     ,
                                            Altitudeh_d  ,
                                            Altitude_d   ,
                                            A            ,
                                            point_local_d,
                                            point_num    ,
                                            1            ,
                                            DeepModel    );
                check_h = false;
                cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
                isnan_check_thor<<< 16, NTH >>>(DivM_d, nv, 3*point_num, check_d);
                isnan_check_thor<<< 16, NTH >>>(divg_Mh_d, nv, 3*point_num, check_d);
                cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
                if(check_h){
                   printf("\n\n Error in NAN check after Thor:DivOp_fast!\n");
                   exit(EXIT_FAILURE);
                }
            }

//          Momentum equation.
            cudaDeviceSynchronize();
            // Updates: Mhs_d
            Momentum_Eq <LN,LN>  <<<NB, NT>>>(Mhs_d      ,
                                              pressures_d,
                                              SlowMh_d   ,
                                              grad_d     ,
                                              Altitude_d ,
                                              diffmh_d     ,
                                              A          ,
                                              func_r_d   ,
                                              times      ,
                                              maps_d     ,
                                              nl_region  ,
                                              DeepModel  );
            // Updates: Mhs_d
            Momentum_Eq_Poles<6><<< 2, 1>>>(Mhs_d        ,
                                            pressures_d  ,
                                            SlowMh_d     ,
                                            grad_d       ,
                                            Altitude_d   ,
                                            DivM_d       ,
                                            A            ,
                                            func_r_d     ,
                                            times        ,
                                            point_local_d,
                                            nv           ,
                                            point_num    ,
                                            DeepModel    );
           check_h = false;
           cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
           isnan_check_thor<<< 16, NTH >>>(Mhs_d, nv, 3*point_num, check_d);
           cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
           if(check_h){
              printf("\n\n Error in NAN check after Thor:Momentum_Eq!\n");
              exit(EXIT_FAILURE);
           }
//          Vertical Momentum
            cudaDeviceSynchronize();

            BENCH_POINT_I_SS(*this, current_step, rk, ns, "Momentum_Eq")
            // Updates: Sp_d, Sd_d
            Prepare_Implicit_Vertical <LN,LN>  <<<NB, NT >>>(Mhs_d         ,
                                                             h_d           ,
                                                             div_d         ,
                                                             Slowpressure_d,
                                                             SlowRho_d     ,
                                                             Sp_d          ,
                                                             Sd_d          ,
                                                             Altitude_d    ,
                                                             Cp            ,
                                                             Rd            ,
                                                             A             ,
                                                             maps_d        ,
                                                             nl_region     ,
                                                             DeepModel     );
           check_h = false;
           cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
           isnan_check_thor<<< 16, NTH >>>(Sp_d, nv, point_num, check_d);
           isnan_check_thor<<< 16, NTH >>>(Sd_d, nv, point_num, check_d);
           cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
           if(check_h){
              printf("\n\n Error in NAN check after Thor:PrepImp!\n");
              exit(EXIT_FAILURE);
           }
            cudaDeviceSynchronize();
            // Updates: Sp_d, Sd_d
            Prepare_Implicit_Vertical_Poles<6><<<2, 1>>> (Mhs_d         ,
                                                          h_d           ,
                                                          div_d         ,
                                                          Slowpressure_d,
                                                          SlowRho_d     ,
                                                          Sp_d          ,
                                                          Sd_d          ,
                                                          Altitude_d    ,
                                                          Cp            ,
                                                          Rd            ,
                                                          A             ,
                                                          point_local_d ,
                                                          point_num     ,
                                                          nv            ,
                                                          DeepModel     );

           check_h = false;
           cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
           isnan_check_thor<<< 16, NTH >>>(Sp_d, nv, point_num, check_d);
           isnan_check_thor<<< 16, NTH >>>(Sd_d, nv, point_num, check_d);
           cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
           if(check_h){
              printf("\n\n Error in NAN check after Thor:PrepImpPoles!\n");
              exit(EXIT_FAILURE);
           }
            cudaDeviceSynchronize();

            // Updates: Whs_d, Ws_d
            Vertical_Eq <<< (point_num / num_th_vertical_eq) + 1,
                num_th_vertical_eq,
                2*num_th_vertical_eq*nvi*sizeof(double) >>> (Whs_d      ,
                                                             Ws_d       ,
                                                             pressures_d,
                                                             h_d        ,
                                                             hh_d       ,
                                                             Rhos_d     ,
                                                             gtil_d     ,
                                                             gtilh_d    ,
                                                             Sp_d       ,
                                                             Sd_d       ,
                                                             SlowWh_d   ,
                                                             Cp         ,
                                                             Rd         ,
                                                             times      ,
                                                             Gravit     ,
                                                             Altitude_d ,
                                                             Altitudeh_d,
                                                             A          ,
                                                             NonHydro   ,
                                                             point_num  ,
                                                             nv         ,
                                                             nvi        ,
                                                             DeepModel  );


             check_h = false;
             cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
             isnan_check_thor<<< 16, NTH >>>(Whs_d, nv+1, point_num, check_d);
             isnan_check_thor<<< 16, NTH >>>(Ws_d, nv, point_num, check_d);
             cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
             if(check_h){
                printf("\n\n Error in NAN check after Thor:VertEq!\n");
                exit(EXIT_FAILURE);
             }

            cudaError_t err = cudaGetLastError();


            // Check device query
            if (err != cudaSuccess)
            {
                printf("%s\n", cudaGetErrorString(err));
            }

//          Pressure and density equations.
            cudaDeviceSynchronize();

            BENCH_POINT_I_SS(*this, current_step, rk, ns, "Vertical_Eq")
            // Updates: pressures_d, Rhos_d
            Density_Pressure_Eqs <LN,LN>  <<<NB, NT >>>(pressures_d,
                                                        pressurek_d,
                                                        Rhos_d     ,
                                                        Rhok_d     ,
                                                        Mhs_d      ,
                                                        Mhk_d      ,
                                                        Whs_d      ,
                                                        Whk_d      ,
                                                        pt_d       ,
                                                        pth_d      ,
                                                        SlowRho_d  ,
                                                        diffpr_d   ,
                                                        div_d      ,
                                                        Altitude_d ,
                                                        Altitudeh_d,
                                                        Cp         ,
                                                        Rd         ,
                                                        A          ,
                                                        P_Ref      ,
                                                        times      ,
                                                        maps_d     ,
                                                        nl_region  ,
                                                        DeepModel  );
           check_h = false;
           cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
           isnan_check_thor<<< 16, NTH >>>(pressures_d, nv, point_num, check_d);
           isnan_check_thor<<< 16, NTH >>>(Rhos_d, nv, point_num, check_d);
           cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
           if(check_h){
              printf("\n\n Error in NAN check after Thor:DensityP!\n");
              exit(EXIT_FAILURE);
           }
            cudaDeviceSynchronize();
            // Updates: pressures_d, Rhos_d
            Density_Pressure_Eqs_Poles<6><<< 2, 1>>>(pressures_d  ,
                                                     pressurek_d  ,
                                                     Rhos_d       ,
                                                     Rhok_d       ,
                                                     Mhs_d        ,
                                                     Mhk_d        ,
                                                     Whs_d        ,
                                                     Whk_d        ,
                                                     pt_d         ,
                                                     pth_d        ,
                                                     SlowRho_d    ,
                                                     diffpr_d     ,
                                                     div_d        ,
                                                     Altitude_d   ,
                                                     Altitudeh_d  ,
                                                     Cp           ,
                                                     Rd           ,
                                                     A            ,
                                                     P_Ref        ,
                                                     times        ,
                                                     point_local_d,
                                                     point_num    ,
                                                     nv           ,
                                                     DeepModel    );
             check_h = false;
             cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
             isnan_check_thor<<< 16, NTH >>>(pressures_d, nv, point_num, check_d);
             isnan_check_thor<<< 16, NTH >>>(Rhos_d, nv, point_num, check_d);
             cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
             if(check_h){
                printf("\n\n Error in NAN check after Thor:DensityPPoles!\n");
                exit(EXIT_FAILURE);
             }

        }

//      Update quantities for the long loop.
        cudaDeviceSynchronize();
        // Updates: Mhk_d, Whk_d, Wk_d, Rhok_d, pressurek_d
        UpdateRK2<<< (point_num / NTH) + 1, NTH >>> (Mhs_d      ,
                                                     Mhk_d      ,
                                                     Whs_d      ,
                                                     Whk_d      ,
                                                     Wk_d       ,
                                                     Rhos_d     ,
                                                     Rhok_d     ,
                                                     pressures_d,
                                                     pressurek_d,
                                                     func_r_d   ,
                                                     Altitude_d ,
                                                     Altitudeh_d,
                                                     point_num  ,
                                                     nv         );

        check_h = false;
        cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
        isnan_check_thor<<< 16, NTH >>>(Mhk_d, nv, 3*point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(Whk_d, nvi, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(Wk_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(Rhok_d, nv, point_num, check_d);
        isnan_check_thor<<< 16, NTH >>>(pressurek_d, nv, point_num, check_d);
        cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
        if(check_h){
           printf("\n\n Error in NAN check after Thor:UpdateRK2!\n");
           exit(EXIT_FAILURE);
        }
        BENCH_POINT_I_S(*this, current_step, rk, "RK2")
    }
//  Update diagnostic variables.
    cudaDeviceSynchronize();

    BENCH_POINT_I(*this, current_step, "END")

    cudaMemcpy(Mh_d       , Mhk_d       , point_num * nv * 3 * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Wh_d       , Whk_d       , point_num * nvi*     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(W_d        , Wk_d        , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(Rho_d      , Rhok_d      , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(pressure_d , pressurek_d , point_num * nv *     sizeof(double), cudaMemcpyDeviceToDevice);
}
//END OF THOR!
