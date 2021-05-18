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

#pragma once

#include "phy_module_base.h"

class radiative_transfer : public phy_module_base
{
public:
    radiative_transfer();
    ~radiative_transfer();

    bool initialise_memory(const ESP &esp, device_RK_array_manager &phy_modules_core_arrays);
    bool initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s);

    bool phy_loop(ESP &                  esp,
                  const SimulationSetup &sim,
                  kernel_diagnostics &   diag,
                  int                    nstep, // Step number
                  double                 time_step);            // Time-step [s]

    bool store(const ESP &esp, storage &s);

    bool store_init(storage &s);

    bool configure(config_file &config_reader);

    virtual bool free_memory();

    void print_config();

    void set_qheat_scaling(const double &scaling) {
        Qheat_scaling = scaling;
    };

    // DEBUG: hack, for debug printout in Alf.
    double *get_debug_qheat_device_ptr() {
        return qheat_d;
    };

    double OLR_tot;
    double ASR_tot;

private:
    // Scaling of Qheat, for slow ramp up or ramp down.
    double Qheat_scaling = 1.0;

    // Config options
    double Tstar_config            = 4520;  // Star effective temperature [K]
    double planet_star_dist_config = 0.015; // Planet-star distance [au]
    double radius_star_config      = 0.667; // Star radius [Rsun]
    double diff_ang_config         = 0.5;   // Diffusivity angle (1 / diffusivity factor): 0.5-1.0
    // double Tint_config       = 0;    // temperature of upward flux coming from the planet's interior
    double albedo_config = 0.18; // Bond albedo
    // double tausw_config      = 532.0;  // Absorption coefficient for the shortwaves
    // double taulw_config      = 1064.0; // Absorption coefficient for the longwaves
    bool   latf_lw_config       = false; // use sin^2(lat) dependence for lw opacity
    double kappa_lw_pole_config = 0.002; // Absorption coefficient for the longwave (poles)
    double n_lw_config          = 2.0;   // power law dependence for unmixed absorbers in LW
    double n_sw_config          = 1.0;   // power law dependence for mixed/unmixed absorbers in SW
    // double f_lw_config       = 0.5;    // fraction of taulw in well-mixed absorber

    // double Csurf_config    = 1e7;   // heat capacity of surface (J K^-1 m^-2)
    bool rt1Dmode_config = false; // 1D mode=all columns are irradiated identically

    int spinup_start_step = -1;
    int spinup_stop_step  = -1;

    int spindown_start_step = -1;
    int spindown_stop_step  = -1;
    // Rad trans
    double Tstar            = 4520;  // Star effective temperature [K]
    double planet_star_dist = 0.015; // Planet-star distance [au]
    double radius_star      = 0.667; // Star radius [Rsun]
    double diff_ang         = 0.5;   // Diffusivity factor: 0.5-1.0
    // double Tint   = 0; // Lower boundary temperature: upward flux coming from the planet's interior
    double albedo   = 0.18;  // Bond albedo
    double kappa_sw = 0.001; // Absorption coefficient for the shortwaves
    double kappa_lw = 0.002; // Absorption coefficient for the longwaves
    double kappa_lw_pole;
    // double kappa_lw_pole = 0.002;
    bool   latf_lw = false;
    double n_lw    = 2.0; // power law dependence for unmixed absorbers in LW
    double n_sw    = 1.0; // power law dependence for mixed/unmixed absorbers in SW
    // double f_lw       = 0.5; // fraction of taulw in well-mixed absorber

    bool rt1Dmode;
    // double  Csurf;
    double *surf_flux_d;


    double incflx;
    //  Arrays used in RT code
    double *fsw_up_d;
    double *fsw_dn_d;
    double *flw_up_d;
    double *flw_dn_d;
    double *tau_d;

    double *tau_h;
    double *fsw_up_h;
    double *fsw_dn_h;
    double *flw_up_h;
    double *flw_dn_h;

    double *insol_h;
    double *insol_d;
    double *insol_ann_h;
    double *insol_ann_d;

    double *qheat_d;
    double *qheat_h;

    double *ASR_d;
    double *OLR_d;


    // picket-fence scheme

    double* gam_P;
    double* gam_V__h;
    double* Beta_V__h;
    double* Beta__h;
    double* gam_1__h;
    double* gam_2__h;
    double* k_IR_2__h;
    double* k_V_3__h;
    double* net_F_h;
    double* AB__h;
    double* Teff;
    
    double  metalicity = 0;
    double* k_IR_2_nv_d;
    double* k_V_3_nv_d;
    double* gam_V_3_d;
    double* gam_1_d;
    double* gam_2_d;
    double* Beta_V_3_d;
    double* Beta_2_d;
    double* net_F_nvi_d;    
    double* AB_d;

    double* lw_net__h ;
    double* sw_net__h ;
    double* dtau__h ;

    //Kitzman working variables                          
    double* tau_Ve__df_e; 
    double* tau_IRe__df_e;
    double* Te__df_e;
    double* be__df_e; 
    double* sw_down__df_e;
    double* sw_down_b__df_e;
    double* sw_up__df_e;
    double* lw_down__df_e;
    double* lw_down_b__df_e;
    double* lw_up__df_e;
    double* lw_up_b__df_e;
    double* lw_net__df_e;
    double* sw_net__df_e;

    // lw_grey_updown_linear working variables
    double* dtau__dff_l;
    double* del__dff_l;
    double* edel__dff_l;
    double* e0i__dff_l;
    double* e1i__dff_l;
    double* Am__dff_l;
    double* Bm__dff_l;
    double* lw_up_g__dff_e;
    double* lw_down_g__dff_e;


    //  These arrays are for temporary usage in RT code
    double *dtemp;
    double *phtemp;
    double *ttemp;
    double *thtemp;
    void    RTSetup(double Tstar_,
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
                    bool   rt1Dmode,
                    double Tmean);
};

void RTSetup_picket_fence(double Tstar_,
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
    double f_lw,
    bool   rt1Dmode_,
    double Tmean);
