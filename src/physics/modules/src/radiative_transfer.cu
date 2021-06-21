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
#include "reduction_add.h"
#include "debug_helpers.h"

#include <stdio.h>
#include <stdlib.h>

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

    double picket_fence_mod = true;

    cudaMalloc((void **)&ASR_d, esp.point_num * sizeof(double));
    cudaMalloc((void **)&OLR_d, esp.point_num * sizeof(double));
    
    

    if (picket_fence_mod){
        //  Rad Transfer
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
        
        cudaMalloc((void **)&phtemp, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&dtemp, esp.nv * esp.point_num * sizeof(double));

        cudaMalloc((void **)&qheat_d, esp.nv * esp.point_num * sizeof(double));
        qheat_h = (double *)malloc(esp.point_num * esp.nv * sizeof(double));

        insol_h = (double *)malloc(esp.point_num * sizeof(double));
        cudaMalloc((void **)&insol_d, esp.point_num * sizeof(double));
       
        
        cudaMalloc((void **)&surf_flux_d, esp.point_num * sizeof(double));

       

        // picket fence parameters

        cudaMalloc((void **)&k_IR_2_nv_d, 2*esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&k_V_3_nv_d, 3*esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&gam_V_3_d, 3 * esp.point_num * sizeof(double));
        cudaMalloc((void **)&gam_1_d,  esp.point_num * sizeof(double));
        cudaMalloc((void **)&gam_2_d,  esp.point_num * sizeof(double));
        cudaMalloc((void **)&Beta_V_3_d, 3* esp.point_num * sizeof(double));
        cudaMalloc((void **)&Beta_2_d, 2 * esp.point_num * sizeof(double));
        cudaMalloc((void **)&net_F_nvi_d, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&AB_d,  esp.point_num * sizeof(double));

        k_IR_2__h = (double *)malloc(2*esp.nv * esp.point_num * sizeof(double));
        k_V_3__h = (double *)malloc(3*esp.nv * esp.point_num * sizeof(double));
        gam_V__h = (double *)malloc(3 * esp.point_num * sizeof(double));
        gam_1__h = (double *)malloc( esp.point_num * sizeof(double));
        gam_2__h = (double *)malloc( esp.point_num * sizeof(double));
        Beta_V__h = (double *)malloc(3* esp.point_num * sizeof(double));
        Beta__h = (double *)malloc(2 * esp.point_num * sizeof(double));
        net_F_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        AB__h = (double *)malloc( esp.point_num * sizeof(double));
        gam_P = (double *)malloc(esp.point_num * sizeof(double));
        Teff = (double *)malloc(esp.point_num * sizeof(double));

        // picket fence parameters     //Kitzman working variables
        cudaMalloc((void **)&tau_Ve__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&tau_IRe__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&Te__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&be__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&sw_down__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&sw_down_b__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&sw_up__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_down__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_down_b__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_up__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_up_b__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_net__df_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&sw_net__df_e, esp.nvi * esp.point_num * sizeof(double));

        lw_net__h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        sw_net__h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
                

        // picket fence parameters     // lw_grey_updown_linear working variables
        cudaMalloc((void **)&dtau__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&del__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&edel__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&e0i__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&e1i__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&Bm__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&Am__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_up_g__dff_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&lw_down_g__dff_e, esp.nvi * esp.point_num * sizeof(double));
        cudaMalloc((void **)&Gp__dff_l, esp.nv * esp.point_num * sizeof(double));
        cudaMalloc((void **)&Bp__dff_l, esp.nv * esp.point_num * sizeof(double));
        


        tau_h    = (double *)malloc(2*esp.nv * esp.point_num * sizeof(double));
        flw_up_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        flw_dn_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        fsw_up_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        fsw_dn_h = (double *)malloc(esp.nvi * esp.point_num * sizeof(double));
        
        
       

    } else {
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

       

    }
    

    return true;
}


bool radiative_transfer::free_memory() {
    double picket_fence_mod = true;

    if (picket_fence_mod) {
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
        
        cudaFree(phtemp);
        cudaFree(dtemp);

        cudaFree(qheat_d);
        free(qheat_h);

        cudaFree(insol_d);
        


        cudaFree(surf_flux_d);
        cudaFree(ASR_d);
        cudaFree(OLR_d);

        // picket fence parameters

        cudaFree(k_IR_2_nv_d);
        cudaFree(k_V_3_nv_d);
        cudaFree(gam_V_3_d);
        cudaFree(gam_1_d);
        cudaFree(gam_2_d);
        cudaFree(Beta_V_3_d);
        cudaFree(Beta_2_d);
        cudaFree(net_F_nvi_d);
        cudaFree(AB_d);
        
        free(k_IR_2__h);
        free(k_V_3__h);
        free(gam_V__h);
        free(gam_1__h);
        free(gam_2__h);
        free(Beta_V__h);
        free(Beta__h);
        free(net_F_h);
        free(AB__h);
        free(gam_P);
        free(Teff);

        // picket fence parameters     //Kitzman working variables

        
        cudaFree(tau_Ve__df_e);        
        cudaFree(tau_IRe__df_e);        
        cudaFree(Te__df_e);        
        cudaFree(be__df_e);        
        cudaFree(sw_down__df_e);        
        cudaFree(sw_down_b__df_e);        
        cudaFree(sw_up__df_e);        
        cudaFree(lw_down__df_e);        
        cudaFree(lw_down_b__df_e);        
        cudaFree(lw_up__df_e);        
        cudaFree(lw_up_b__df_e);        
        cudaFree(lw_net__df_e);        
        cudaFree(sw_net__df_e);

        free(lw_net__h);
        free(sw_net__h);
                        

        // picket fence parameters     // lw_grey_updown_linear working variables
        cudaFree(dtau__dff_l);
        cudaFree(del__dff_l);
        cudaFree(edel__dff_l);
        cudaFree(e0i__dff_l);
        cudaFree(e1i__dff_l);
        cudaFree(Bm__dff_l);
        cudaFree(Am__dff_l);
        cudaFree(lw_up_g__dff_e);
        cudaFree(lw_down_g__dff_e);
        cudaFree(Gp__dff_l);
        cudaFree(Bp__dff_l);
        

        free(tau_h);
        
        

    } else {
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
        cudaFree(ASR_d);
        cudaFree(OLR_d);

    }
    

    return true;
}

bool radiative_transfer::initial_conditions(const ESP &            esp,
                                            const SimulationSetup &sim,
                                            storage *              s) {
                                            
    cuda_check_status_or_exit(__FILE__, __LINE__);

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

    double picket_fence_mod = true;

    if (picket_fence_mod) {
    
        
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
            //esp.f_lw,
            rt1Dmode_config,
            sim.Tmean);
            
            
    } else {
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
            //esp.f_lw,
            rt1Dmode_config,
            sim.Tmean);

        cudaMemset(surf_flux_d, 0, sizeof(double) * esp.point_num);

    }

    

    

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

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


// Calculates the Bond Albedo according to Parmentier et al. (2015) expression
void Bond_Parmentier(double Teff0, double grav, double& AB) {
    // dependcies
    //// pow from math
    //// log10 from math        


    // Input:
    // Teff0 - Atmospheric profile effective temperature [K] with zero albedo
    // grav - Surface gravity of planet [m s-2]

    // Call by reference (Input&Output):
    // AB - Bond albedo

    // work variables
    double a = 0.0, b = 0.0;

    // start operations

    if (Teff0 <= 250.0)
    {
        a = ((double)-0.335) * pow(grav, ((double)0.070));
        b = 0.0;
    }
    else if (Teff0 > 250.0 && Teff0 <= 750.0)
    {
        a = -0.335 * pow(grav, ((double)0.070)) + 2.149 * pow(grav, ((double)0.135));
        b = -0.896 * pow(grav, ((double)0.135));
    }
    else if (Teff0 > 750.0 && Teff0 < 1250.0)
    {
        a = -0.335 * pow(grav, ((double)0.070)) - 0.428 * pow(grav, ((double)0.135));
        b = 0.0;
    }
    else if (Teff0 >= 1250.0)
    {
        a = 16.947 - ((double)3.174) * pow(grav, ((double)0.070)) - 4.051 *
            pow(grav, ((double)0.135));
        b = -5.472 + ((double)0.917) * pow(grav, ((double)0.070)) + 1.170 *
            pow(grav, ((double)0.135));
    }

    // Final Bond Albedo expression
    AB = pow(((double)10.0), (a + b * log10(Teff0)));

}


///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Calculates 3 band grey visual gamma values and 2 picket fence IR gamma values
    // according to the coefficents and equations in:
    // Parmentier & Menou (2014) and Parmentier et al. (2015)
    // NOTE: This does not calculate the opacity - call k_Ross_Freedman for that
    void gam_Parmentier(int nCol, int nLev, double *Teff, int table_num, double *gam_V, double *Beta_V,
    double *Beta, double *gam_1, double *gam_2 , double *gam_P) {
    // dependcies
    //// pow from math
    //// log10 from math        


    // Input:
    // Teff - Effective temperature [K] (See Parmentier papers for various ways to calculate this)
    // for non-irradiated atmosphere Teff = Tint
    // table_num - Table selection from Parmentier et al. (2015): 1 = w. TiO/VO, 2 = w.o. TiO/VO

    // Call by reference (Input&Output):
    // gam_V(3) - gamma ratio for 3 visual bands (gam_V = kV_Ross/kIR_Ross)
    // beta_V(3) - fraction of total incident stellar flux in band (1/3 for Parmentier values)
    // Beta - equilvalent bandwidth for picket fence IR model
    // gam_1 - gamma ratio for IR band 1 (gam_1 = kIR_1/kIR_Ross)
    // gam_2 - gamma ratio for IR band 2 (gam_2 = kIR_2/kIR_Ross)
    // gam_P - gamma ratio for Planck mean (gam_P = kIR_Planck/kIR_Ross)
    // tau_lim - tau limit variable (usually for IC system)

    // work variables
    double  R = 0;
    double aP = 0;
    double bP = 0;
    double cP = 0;
    double aV1 = 0, bV1 = 0, aV2 = 0, bV2 = 0, aV3 = 0, bV3 = 0;
    double aB = 0, bB = 0;
    double l10T = 0, l10T2 = 0, RT = 0;
    int i;

   

    // start operations

    for (int id = 0; id < nCol; id++) {

        // Log 10 T_eff variables
        l10T = log10(Teff[id]);
        l10T2 = pow(l10T, 2.0);

        if (table_num == 1) {
            // First table in Parmentier et al. (2015) w. TiO/VO
            // Start large if statements with visual band and Beta coefficents
            if (Teff[id] <= 200.0)
            {
                aV1 = -5.51; bV1 = 2.48;
                aV2 = -7.37; bV2 = 2.53;
                aV3 = -3.03; bV3 = -0.20;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 200.0 && Teff[id] <= 300.0)
            {
                aV1 = 1.23; bV1 = -0.45;
                aV2 = 13.99; bV2 = -6.75;
                aV3 = -13.87; bV3 = 4.51;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 300.0 && Teff[id] <= 600.0)
            {
                aV1 = 8.65; bV1 = -3.45;
                aV2 = -15.18; bV2 = 5.02;
                aV3 = -11.95; bV3 = 3.74;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 600.0 && Teff[id] <= 1400.0)
            {
                aV1 = -12.96; bV1 = 4.33;
                aV2 = -10.41; bV2 = 3.31;
                aV3 = -6.97; bV3 = 1.94;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 1400.0 && Teff[id] < 2000.0)
            {
                aV1 = -23.75; bV1 = 7.76;
                aV2 = -19.95; bV2 = 6.34;
                aV3 = -3.65; bV3 = 0.89;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] >= 2000.0)
            {
                aV1 = 12.65; bV1 = -3.27;
                aV2 = 13.56; bV2 = -3.81;
                aV3 = -6.02; bV3 = 1.61;
                aB = 6.21; bB = -1.63;
            }

            // gam_P coefficents
            aP = -2.36;
            bP = 13.92;
            cP = -19.38;
        }
        else if (table_num == 2)
        {
            // ! Appendix table from Parmentier et al. (2015) - without TiO and VO
            if (Teff[id] <= 200.0)
            {
                aV1 = -5.51; bV1 = 2.48;
                aV2 = -7.37; bV2 = 2.53;
                aV3 = -3.03; bV3 = -0.20;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 200.0 && Teff[id] <= 300.0)
            {
                aV1 = 1.23; bV1 = -0.45;
                aV2 = 13.99; bV2 = -6.75;
                aV3 = -13.87; bV3 = 4.51;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 300.0 && Teff[id] <= 600.0)
            {
                aV1 = 8.65; bV1 = -3.45;
                aV2 = -15.18; bV2 = 5.02;
                aV3 = -11.95; bV3 = 3.74;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 600.0 && Teff[id] <= 1400.0)
            {
                aV1 = -12.96; bV1 = 4.33;
                aV2 = -10.41; bV2 = 3.31;
                aV3 = -6.97; bV3 = 1.94;
                aB = 0.84; bB = 0.0;
            }
            else if (Teff[id] > 1400.0 && Teff[id] < 2000.0)
            {
                aV1 = -1.68; bV1 = 0.75;
                aV2 = 6.96; bV2 = -2.21;
                aV3 = 0.02; bV3 = -0.28;
                aB = 3.0; bB = -0.69;
            }
            else if (Teff[id] >= 2000.0)
            {
                aV1 = 10.37; bV1 = -2.91;
                aV2 = -2.4; bV2 = 0.62;
                aV3 = -16.54; bV3 = 4.74;
                aB = 3.0; bB = -0.69;
            }

            // gam_P coefficents
            if (Teff[id] <= 1400.0)
            {
                aP = -2.36;
                bP = 13.92;
                cP = -19.38;
            }
            else
            {
                aP = -12.45;
                bP = 82.25;
                cP = -134.42;
            }
        }

        // Calculation of all values
        // Visual band gamma
        gam_V[id*3 + 0] = pow(((double)10.0), (aV1 + bV1 * l10T));
        gam_V[id*3 + 1] = pow(((double)10.0), (aV2 + bV2 * l10T));
        gam_V[id*3 + 2] = pow(((double)10.0), (aV3 + bV3 * l10T));



        // Visual band fractions
        for (i = 0; i < 3; i++)
        {
            Beta_V[id*3 + i] = ((double)1.0) / ((double)3.0);
        }

        // gamma_Planck - if < 1 then make it grey approximation (k_Planck = k_Ross, gam_P = 1)
        gam_P[id] = pow(((double)10.0), (aP * l10T2 + bP * l10T + cP));
        if (gam_P[id] < 1.0000001)
        {
            gam_P[id] = 1.0000001;
        }

        // equivalent bandwidth value
        Beta[id*2 + 0] = aB + bB * l10T;
        Beta[id*2 + 1] = (1.0) - Beta[id*2 + 0];

        // IR band kappa1/kappa2 ratio - Eq. 96 from Parmentier & Menou (2014)
        RT = (gam_P[id] - 1.0) / (2.0 * Beta[id*2 + 0] * Beta[id*2 + 1]);
        R = 1.0 + RT + sqrt(pow(RT, 2.0) + RT);

        // gam_1 and gam_2 values - Eq. 92, 93 from Parmentier & Menou (2014)
        gam_1[id] = Beta[id*2 + 0] + R - Beta[id*2 + 0] * R;
        gam_2[id] = gam_1[id] / R;
    }

    

    

}

bool radiative_transfer::phy_loop(ESP &                  esp,
                                  const SimulationSetup &sim,
                                  kernel_diagnostics &   diag,
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
        double picket_fence_mod = true;
        double Tirr;
        //
        //  Number of threads per block.
        const int NTH = 256;

        //  Specify the block sizes.
        dim3 NB((esp.point_num / NTH) + 1, esp.nv, 1);
        dim3 NBRT((esp.point_num / NTH) + 1, 1, 1);

        if (picket_fence_mod){
                
            for (int c = 0; c <  esp.point_num; c++){
                // Parmentier opacity profile parameters - first get Bond albedo
                Tirr = Tstar*pow((radius_star/planet_star_dist), 0.5);
                
                Teff[c] = pow( (pow(esp.Tint, 4.0) +
                    (1.0 / sqrt(3.0)) *
                    pow(Tirr, 4.0) ), 0.25);

                Bond_Parmentier(Teff[c], sim.Gravit, AB__h[c]);
                
               
                // Recalculate Teff and then find parameters
                if (esp.insolation.get_host_cos_zenith_angles()[c] >= 0.0)
                {
                    Teff[c] = pow( (pow(esp.Tint, 4.0) +
                        (1.0 - AB__h[c]) * esp.insolation.get_host_cos_zenith_angles()[c] * pow(Tirr, 4.0) ), 0.25);
                } else {
                    Teff[c] = pow( pow(esp.Tint, 4.0) + 0.0, 0.25);
                }
                
            }
            
            
            
            gam_Parmentier(esp.point_num,
                esp.nv,
                Teff,
                2,
                gam_V__h,
                Beta_V__h,
                Beta__h,
                gam_1__h,
                gam_2__h,
                gam_P);
                
            
                
            bool cudaStatus;
            cudaStatus = cudaMemcpy(gam_V_3_d, gam_V__h, 3*esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "gam_V_3_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }
            cudaStatus = cudaMemcpy(Beta_V_3_d, Beta_V__h, 3*esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "Beta_V_3_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }
            cudaStatus = cudaMemcpy(Beta_2_d, Beta__h, 2*esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "Beta_2_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }
            cudaStatus = cudaMemcpy(gam_1_d, gam_1__h, esp.point_num *  sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "gam_1_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }
            cudaStatus = cudaMemcpy(gam_2_d, gam_2__h, esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "gam_2_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }
            cudaStatus = cudaMemcpy(AB_d, AB__h, esp.point_num * sizeof(double), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "AB_d cudaMemcpyHostToDevice failed!");
                //goto Error;
            }

            // check for error
            cudaError_t error = cudaGetLastError();
            if(error != cudaSuccess)
            {
                // print the CUDA error message and exit
                printf("CUDA error: %s\n", cudaGetErrorString(error));
                exit(-1);
            }            
            
                      

            rtm_picket_fence<<<NBRT, NTH>>>(esp.pressure_d,
                esp.temperature_d,
                esp.Rho_d,
                sim.Gravit,
                esp.Cp_d,
                esp.lonlat_d,
                esp.Altitude_d,
                esp.Altitudeh_d,
                esp.insolation.get_r_orb(),
                phtemp,
                dtemp,
                time_step,
                esp.Tint,
                albedo,
                kappa_lw,
                latf_lw,
                kappa_lw_pole,
                incflx,
                esp.point_num,
                esp.nv,
                esp.nvi,
                sim.A,
                esp.insolation.get_device_cos_zenith_angles(),
                insol_d,
                esp.surface,
                esp.Tsurface_d,
                surf_flux_d,
                esp.areasT_d,
                ASR_d,
                OLR_d,
                esp.profx_Qheat_d,
                qheat_d,
                esp.Rd_d,
                Qheat_scaling,
                metalicity,
                k_IR_2_nv_d,
                k_V_3_nv_d,
                gam_V_3_d,
                gam_1_d,
                gam_2_d,
                Beta_V_3_d,
                Beta_2_d,
                net_F_nvi_d,    
                AB_d,
                tau_Ve__df_e, //Kitzman working variables
                tau_IRe__df_e, 
                Te__df_e, 
                be__df_e, 
                sw_down__df_e,
                sw_down_b__df_e,
                sw_up__df_e,
                lw_down__df_e, 
                lw_down_b__df_e,
                lw_up__df_e,
                lw_up_b__df_e,
                lw_net__df_e,
                sw_net__df_e,
                dtau__dff_l, // lw_grey_updown_linear working variables
                del__dff_l, 
                edel__dff_l,
                e0i__dff_l, 
                e1i__dff_l,
                Am__dff_l, 
                Bm__dff_l,
                lw_up_g__dff_e, 
                lw_down_g__dff_e,
                Gp__dff_l,
                Bp__dff_l,

                rt1Dmode,
                sim.DeepModel                                        

                );



                

        } else {

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
                kappa_sw,
                kappa_lw,
                latf_lw,
                kappa_lw_pole,
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
                esp.dTsurf_dt_d,
                surf_flux_d,
                esp.areasT_d,
                ASR_d,
                OLR_d,
                esp.profx_Qheat_d,
                qheat_d,
                esp.Rd_d,
                Qheat_scaling,
                sim.gcm_off,
                rt1Dmode,
                sim.DeepModel);

               

        }

        
        

        
        ASR_tot = gpu_sum_on_device<1024>(ASR_d, esp.point_num);
        OLR_tot = gpu_sum_on_device<1024>(OLR_d, esp.point_num);
        
        cuda_check_status_or_exit(__FILE__, __LINE__);

        if (nstep * time_step < (2 * M_PI / esp.insolation.get_mean_motion())) {
            // stationary orbit/obliquity
            // calculate annually average of insolation for the first orbit

            double picket_fence_mod = true;

            if (picket_fence_mod) {

            } else {
                annual_insol<<<NBRT, NTH>>>(insol_ann_d, insol_d, nstep, esp.point_num);
            }
            
        }
    }
    cudaDeviceSynchronize();
    return true;
}

bool radiative_transfer::configure(config_file &config_reader) {
    double picket_fence_mod = true;

    if (picket_fence_mod) {
    
        cuda_check_status_or_exit(__FILE__, __LINE__);
        
        // basic star-planet properties
        config_reader.append_config_var("Tstar", Tstar_config, Tstar_config);
        config_reader.append_config_var(
            "planet_star_dist", planet_star_dist_config, planet_star_dist_config);
        config_reader.append_config_var("radius_star", radius_star_config, radius_star_config);
        //config_reader.append_config_var("diff_ang", diff_ang_config, diff_ang_config);
        // config_reader.append_config_var("Tint", Tint_config, Tint_config);
        config_reader.append_config_var("albedo", albedo_config, albedo_config);
        // config_reader.append_config_var("tausw", tausw_config, tausw_config);
        // config_reader.append_config_var("taulw", taulw_config, taulw_config);

        // options for latitude dependence in longwave opacity
        config_reader.append_config_var("latf_lw", latf_lw_config, latf_lw_config);
        config_reader.append_config_var("kappa_lw_pole", kappa_lw_pole_config, kappa_lw_pole_config);

       
        // config_reader.append_config_var("f_lw", f_lw_config, f_lw_config);


        config_reader.append_config_var("rt1Dmode", rt1Dmode_config, rt1Dmode_config);

        // spin up spin down
        config_reader.append_config_var("dgrt_spinup_start", spinup_start_step, spinup_start_step);
        config_reader.append_config_var("dgrt_spinup_stop", spinup_stop_step, spinup_stop_step);
        config_reader.append_config_var(
            "dgrt_spindown_start", spindown_start_step, spindown_start_step);
        config_reader.append_config_var("dgrt_spindown_stop", spindown_stop_step, spindown_stop_step);
        
        cuda_check_status_or_exit(__FILE__, __LINE__);


    } else {
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
    }


    

    return true;
}

bool radiative_transfer::store(const ESP &esp, storage &s) {

    double picket_fence_mod = true;

    if (picket_fence_mod) {
        cudaMemcpy(insol_h, insol_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(insol_h, esp.point_num, "/insol", "W m^-2", "insolation (instantaneous)");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
    
        cudaMemcpy(
            lw_net__h, lw_net__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(lw_net__h, esp.nvi * esp.point_num, "/lw_net__h", "W m^-2", "net long-wave flux (LW)");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
    
        cudaMemcpy(
            sw_net__h, sw_net__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(sw_net__h, esp.nvi * esp.point_num, "/sw_net__h", "W m^-2", "net short-wave flux (SW)");

        cudaMemcpy(
            flw_up_h, lw_up__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(flw_up_h, esp.nvi * esp.point_num, "/flw_up", "W m^-2", "upward flux (LW)");
    
        cudaMemcpy(
            fsw_up_h, sw_up__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(fsw_up_h, esp.nvi * esp.point_num, "/fsw_up", "W m^-2", "upward flux (SW)");
    
        cudaMemcpy(
            flw_dn_h, lw_down__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(flw_dn_h, esp.nvi * esp.point_num, "/flw_dn", "W m^-2", "downward flux (LW)");
    
        cudaMemcpy(
            fsw_dn_h, sw_down__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(fsw_dn_h, esp.nvi * esp.point_num, "/fsw_dn", "W m^-2", "downward flux (SW)");

        cuda_check_status_or_exit(__FILE__, __LINE__); 

        //cudaMemcpy(tau_h, tau_Ve__df_e, esp.nvi * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(tau_h,
                       2*esp.nv * esp.point_num,
                       "/tau",
                       " ",
                       "optical depth across each layer (not total optical depth)");
    
    
        cuda_check_status_or_exit(__FILE__, __LINE__);             
    
        cudaMemcpy(qheat_h, qheat_d, esp.nv * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(qheat_h, esp.nv * esp.point_num, "/DGQheat", " ", "Double Gray Qheat");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
    
        // cudaMemcpy(Tsurface_h, Tsurface_d, esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        // s.append_table(Tsurface_h, esp.point_num, "/Tsurface", "K", "surface temperature");
    
        s.append_value(Qheat_scaling, "/dgrt_qheat_scaling", " ", "Qheat scaling applied to DG");
    
        s.append_value(ASR_tot, "/ASR", "W", "Absorbed Shortwave Radiation (global total)");
        s.append_value(OLR_tot, "/OLR", "W", "Outgoing Longwave Radiation (global total)");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);

        cudaMemcpy(k_IR_2__h, k_IR_2_nv_d, 2*esp.nv * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(k_IR_2__h, 2*esp.nv * esp.point_num, "/k_IR_2__h", " ", "kappa for two IR bands");
        cudaMemcpy(k_V_3__h, k_V_3_nv_d,3* esp.nv * esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
        s.append_table(k_V_3__h, 3*esp.nv * esp.point_num, "/k_V_3__h", " ", "kappa for three V bands");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);

              


    } else {
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
    
        s.append_value(ASR_tot, "/ASR", "W", "Absorbed Shortwave Radiation (global total)");
        s.append_value(OLR_tot, "/OLR", "W", "Outgoing Longwave Radiation (global total)");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);
    }
   

    return true;
}

bool radiative_transfer::store_init(storage &s) {
    double picket_fence_mod = true;

    if (picket_fence_mod) {
        s.append_value(Tstar, "/Tstar", "K", "Temperature of host star");
        // s.append_value(Tint, "/Tint", "K", "Temperature of interior heat flux");
        s.append_value(
            planet_star_dist_config, "/planet_star_dist", "au", "distance b/w host star and planet");
        s.append_value(radius_star_config, "/radius_star", "R_sun", "radius of host star");
        
        s.append_value(albedo, "/albedo", "-", "bond albedo of planet");
        //  s.append_value(kappa_sw, "/kappa_sw", "-", "gray opacity of shortwave");
        //  s.append_value(kappa_lw, "/kappa_lw", "-", "gray opacity of longwave");
    
        s.append_value(latf_lw ? 1.0 : 0.0, "/latf_lw", "-", "use lat dependent opacity");
        s.append_value(kappa_lw_pole, "/kappa_lw_pole", "-", "gray opacity of longwave at poles");
       
        // s.append_value(f_lw, "/f_lw", "-", "fraction of taulw in well-mixed absorber");
        
        cuda_check_status_or_exit(__FILE__, __LINE__);

    } else {
        s.append_value(Tstar, "/Tstar", "K", "Temperature of host star");
        // s.append_value(Tint, "/Tint", "K", "Temperature of interior heat flux");
        s.append_value(
            planet_star_dist_config, "/planet_star_dist", "au", "distance b/w host star and planet");
        s.append_value(radius_star_config, "/radius_star", "R_sun", "radius of host star");
        s.append_value(diff_ang, "/diff_ang", "-", "diffusivity factor");
        s.append_value(albedo, "/albedo", "-", "bond albedo of planet");
        //  s.append_value(kappa_sw, "/kappa_sw", "-", "gray opacity of shortwave");
        //  s.append_value(kappa_lw, "/kappa_lw", "-", "gray opacity of longwave");
    
        s.append_value(latf_lw ? 1.0 : 0.0, "/latf_lw", "-", "use lat dependent opacity");
        s.append_value(kappa_lw_pole, "/kappa_lw_pole", "-", "gray opacity of longwave at poles");
        s.append_value(n_lw, "/n_lw", "-", "power law exponent for unmixed absorber in LW");
        s.append_value(n_sw, "/n_sw", "-", "power law exponent for mixed/unmixed absorber in SW");
        // s.append_value(f_lw, "/f_lw", "-", "fraction of taulw in well-mixed absorber");
    }
   


    return true;
}

void radiative_transfer::RTSetup(double Tstar_,
                                 double planet_star_dist_,
                                 double radius_star_,
                                 double diff_ang_,
                                 double P_Ref,
                                 double Gravit,
                                 double albedo_,
                                 double kappa_sw_,
                                 double kappa_lw_,
                                 bool   latf_lw_,
                                 double kappa_lw_pole_,
                                 double n_lw_,
                                 double n_sw_,
                                 //double f_lw,
                                 bool   rt1Dmode_,
                                 double Tmean) {

    double bc = 5.6703744191844314e-08; // Stefan–Boltzmann constant [W m−2 K−4]

    Tstar            = Tstar_;
    planet_star_dist = planet_star_dist_ * 149597870.7; //conv to km
    radius_star      = radius_star_ * 695700;           //conv to km
    diff_ang         = diff_ang_;
    // Tint             = Tint_;
    albedo = albedo_;
    //tausw      = kappa_sw * P_Ref / Gravit;
    //taulw      = kappa_lw * P_Ref / (f_lw * Gravit);
    //taulw_pole = kappa_lw_pole * P_Ref / (f_lw * Gravit);
    kappa_sw      = kappa_sw_;
    kappa_lw      = kappa_lw_;
    kappa_lw_pole = kappa_lw_pole_;

    latf_lw = latf_lw_;
    n_sw    = n_sw_;
    n_lw    = n_lw_;
    //f_lw             = f_lw_;

    double resc_flx = pow(radius_star / planet_star_dist, 2.0);
    incflx          = resc_flx * bc * Tstar * Tstar * Tstar * Tstar;

    rt1Dmode = rt1Dmode_;
    
    cuda_check_status_or_exit(__FILE__, __LINE__);
}

