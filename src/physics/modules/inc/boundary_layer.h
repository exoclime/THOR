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
// Method: Boundary layer (surface friction) physics module
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

#include "diagnostics.h"
#include "phy_module_base.h"

#define bl_type_default "RayleighHS"

#define LARGERiB 1e8

#define KVONKARMAN 0.4

//for non-local scheme, not currently used
#define GAMMA_M 15.0
#define GAMMA_H 15.0

//tuning parameters for thermals (non-local scheme)
#define A_THERM 0.0
#define B_THERM 0.0

//tuning parameters for local diff scheme
// #define TRANSITION_HEIGHT -1
// #define ABL_ASYM_LEN 150.0
// #define FREE_ASYM_LEN 30.0
#define E_MIN_MIX 0 //no idea what this should be!

enum boundary_layer_types { RAYLEIGHHS = 0, LOCALMIXL = 1 };


class boundary_layer : public phy_module_base
{
public:
    boundary_layer();
    ~boundary_layer();

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

private:
    // Config options
    boundary_layer_types bl_type;
    string               bl_type_str;

    // parameters for RayleighHS scheme
    double surf_drag_config = 1.0 / 86400.0; // surface drag coefficient
    double bl_sigma_config  = 0.7;           // sigma coord of boundary layer (% of surf pressure)
    double surf_drag;
    double bl_sigma;

    double asl_transition_height_config =
        -1; //transition height for asymptotic scale length (from BL to free atmos)
            //  (set to -1 for no transition->use abl_asym_len everywhere)
    double abl_asym_len_config  = 150.0; //asymptotic scale length (meters) for BL
    double free_asym_len_config = 30.0;  //asympt. scale length (m) for free atmos

    double asl_transition_height;
    double abl_asym_len;
    double free_asym_len;

    double *cpr_tmp, *dpr_tmp;
    double  zbl; // altitude of transition from BL to free atmosph (ekman scheme)

    int *bl_top_lev_d; // index of highest level (center) inside BL
    int *bl_top_lev_h; // index of highest level (center) inside BL

    double *bl_top_height_d; // height of bl
    double *bl_top_height_h; // height of bl

    double *F_sens_d; //surface to atmos sensible heat flux
    double *F_sens_h;

    double *KM_d; // momentum diffusivity (turbulence)
    double *KM_h;
    double *Rho_int_d; //density at interfaces
    double *p_int_d;   // pressure at interfaces

    double *KH_d; // heat diffusivity (turbulence)
    double *KH_h;

    double *CM_d; // surface drag coeff
    double *CM_h;
    double *CH_d; // surface heat-transfer coeff
    double *CH_h;
    double *vh_lowest_d; //speed of lowest layer
    double *pt_surf_d;   //pt of surface

    double *RiGrad_d; //gradient Ri number
    double *RiGrad_h; //gradient Ri number

    double Ri_crit_config = 1.0; //not used in local scheme right now, but could be
    double z_rough_config = 3.21e-5;

    double Ri_crit; //critical Richardson number (not used currently)
    double z_rough; // roughness length (scale for momentum)

    void BLSetup(const ESP &            esp,
                 const SimulationSetup &sim,
                 int                    bl_type_,
                 double                 surf_drag_,
                 double                 bl_sigma_,
                 double                 Ri_crit_,
                 double                 z_rough_,
                 double                 asl_transition_height_,
                 double                 abl_asym_len_,
                 double                 free_asym_len_);
};

__global__ void rayleighHS(double *Mh_d,
                           double *pressure_d,
                           double *Rho_d,
                           double *Altitude_d,
                           double  surf_drag,
                           double  bl_sigma,
                           double  Gravit,
                           double  time_step,
                           double  A,
                           int     num,
                           bool    GravHeightVar);

__global__ void Momentum_Diff_Impl(double *      Mh_d,
                                   double *      pressure_d,
                                   double *      Rho_d,
                                   double *      Altitude_d,
                                   double *      Altitudeh_d,
                                   double *      cpr_tmp,
                                   double *      dpr_tmp,
                                   double *      KM_d,
                                   double *      Rho_int_d,
                                   double        time_step,
                                   double        A,
                                   int           num,
                                   int           nv,
                                   int *         bl_top_lev_d,
                                   bool          DeepModel,
                                   unsigned int *diagnostics_flag,
                                   diag_data *   diagnostics_data);

__global__ void Heat_Diff_Impl_EnergyEq(double *      pt_d,
                                        double *      pressure_d,
                                        double *      temperature_d,
                                        double *      Rho_d,
                                        double *      Altitude_d,
                                        double *      Altitudeh_d,
                                        double *      Tsurface_d,
                                        double *      cpr_tmp,
                                        double *      dpr_tmp,
                                        double *      KH_d,
                                        double *      Rho_int_d,
                                        double *      pt_surf_d,
                                        double *      p_int_d,
                                        double *      F_sens_d,
                                        double        time_step,
                                        double        Rd,
                                        double        Cp,
                                        double        P_Ref,
                                        double        Csurf,
                                        double        A,
                                        int           num,
                                        int           nv,
                                        int *         bl_top_lev_d,
                                        bool          DeepModel,
                                        unsigned int *diagnostics_flag,
                                        diag_data *   diagnostics_data);

__device__ double stability_fM_u(double CN, double Ri, double z, double z_rough);
__device__ double stability_fM_u(double CN, double Ri, double z, double z_rough);
__device__ double stability_f_s(double Ri);

__global__ void CalcGradRi(double *pressure_d,
                           double *Rho_d,
                           double *Mh_d,
                           double *Tsurface_d,
                           double *pt_d,
                           double *Altitude_d,
                           double *Altitudeh_d,
                           double  Rd,
                           double  Cp,
                           double  P_Ref,
                           double  Gravit,
                           double  Ri_crit,
                           double  z_rough,
                           double *RiGrad_d,
                           double *pt_surf_d,
                           double *p_int_d,
                           double *CM_d,
                           double *CH_d,
                           double *Rho_int_d,
                           double *KM_d,
                           double *KH_d,
                           double *F_sens_d,
                           double  asl_transition_height,
                           double  abl_asym_len,
                           double  free_asym_len,
                           double  A,
                           int     num,
                           int     nv,
                           bool    GravHeightVar);

__global__ void FreeAtmosCutOff(double *KH_d,
                                double *KM_d,
                                double *RiGrad_d,
                                double *Altitudeh_d,
                                double  Ri_crit,
                                int     num,
                                int     nv);
