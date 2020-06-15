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

#define bl_type_default "RayleighHS"

#define LARGERiB 1e8

#define KVONKARMAN 0.4

enum boundary_layer_types { RAYLEIGHHS = 0, MONINOBUKHOV = 1, EKMANSPIRAL = 2 };


class boundary_layer : public phy_module_base
{
public:
    boundary_layer();
    ~boundary_layer();

    bool initialise_memory(const ESP &esp, device_RK_array_manager &phy_modules_core_arrays);
    bool initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s);

    bool phy_loop(ESP &                  esp,
                  const SimulationSetup &sim,
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

    double surf_drag_config = 1.0 / 86400.0; // surface drag coefficient
    double bl_sigma_config  = 0.7;           // sigma coord of boundary layer (% of surf pressure)

    double surf_drag;
    double bl_sigma;
    // double *dvdz_tmp;
    double *d2vdz2_tmp;
    double *atmp, *btmp, *ctmp, *cpr_tmp, *dtmp, *dpr_tmp;
    double  zbl; // altitude of transition from BL to free atmosph (ekman scheme)

    int *   bl_top_lev_d;    // index of highest level inside BL
    int *   bl_top_lev_h;    // index of highest level inside BL
    double *bl_top_height_d; // height of bl
    double *bl_top_height_h; // height of bl

    double *RiB_d;  // bulk Richardson number
    double *RiB_h;  // bulk Richardson number
    double *zeta_d; // m-o stability parameter

    double *KM_d; // momentum diffusivity (turbulence)
    double *KM_h;

    double *KH_d; // heat diffusivity (turbulence)
    double *KH_h;

    double *CD_d; // surface drag coeff
    double *CD_h;
    double *CH_d; // surface heat-transfer coeff
    double *CH_h;
    double *vh_lowest_d; //speed of lowest layer
    double *pt_surf_d;   //pt of surface
    double *p_surf_d;

    double Ri_crit_config      = 1.0;
    double z_rough_config      = 3.21e-5;
    double z_therm_config      = 3.21e-5;
    double f_surf_layer_config = 0.1;

    double Ri_crit;      //critical Richardson number
    double z_rough;      // roughness length (scale for momentum)
    double z_therm;      // thermal "roughness" length
    double f_surf_layer; //fraction of BL in surface layer

    void BLSetup(const ESP &            esp,
                 const SimulationSetup &sim,
                 int                    bl_type_,
                 double                 surf_drag_,
                 double                 bl_sigma_,
                 double                 Ri_crit_,
                 double                 z_rough_,
                 double                 z_therm_,
                 double                 f_surf_layer_);
};

__global__ void rayleighHS(double *Mh_d,
                           double *pressure_d,
                           double *Rho_d,
                           double *Altitude_d,
                           double  surf_drag,
                           double  bl_sigma,
                           double  Gravit,
                           double  time_step,
                           int     num);

__global__ void ConstKMEkman(double *Mh_d,
                             double *pressure_d,
                             double *Rho_d,
                             double *Altitude_d,
                             double *Altitudeh_d,
                             double *d2vdz2_tmp,
                             double  KMconst,
                             double  zbl,
                             double  time_step,
                             int     num,
                             int     nv);

__global__ void Momentum_Diff_Impl(double *Mh_d,
                                   double *pressure_d,
                                   double *Rho_d,
                                   double *Altitude_d,
                                   double *Altitudeh_d,
                                   double *atmp,
                                   double *btmp,
                                   double *ctmp,
                                   double *cpr_tmp,
                                   double *dtmp,
                                   double *dpr_tmp,
                                   double *KM_d,
                                   double  time_step,
                                   int     num,
                                   int     nv,
                                   int *   bl_top_lev_d);

__global__ void Heat_Diff_Impl(double *pt_d,
                               double *pressure_d,
                               double *Rho_d,
                               double *Altitude_d,
                               double *Altitudeh_d,
                               double *Tsurface_d,
                               double *atmp,
                               double *btmp,
                               double *ctmp,
                               double *cpr_tmp,
                               double *dtmp,
                               double *dpr_tmp,
                               double *KH_d,
                               double *pt_surf_d,
                               double *p_surf_d,
                               double  time_step,
                               double  Rd,
                               double  Cp,
                               double  P_Ref,
                               double  Csurf,
                               int     num,
                               int     nv,
                               int *   bl_top_lev_d);

__global__ void CalcRiB(double *pressure_d,
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
                        double  z_therm,
                        double *RiB_d,
                        int *   bl_top_lev_d,
                        double *bl_top_height_d,
                        double *pt_surf_d,
                        double *p_surf_d,
                        double *CD_d,
                        double *CH_d,
                        double *zeta_d,
                        double *vh_lowest_d,
                        int     num,
                        int     nv);

__device__ double RiB_2_zeta(double RiB, double Ri_crit, double z_z0, double z_zT);

__global__ void CalcKM_KH(double *RiB_d,
                          double *zeta_d,
                          double *CD_d,
                          double *CH_d,
                          double *bl_top_height_d,
                          int *   bl_top_lev_d,
                          double *vh_lowest_d,
                          double *Altitude_d,
                          double *Altitudeh_d,
                          double  Ri_crit,
                          double  z_rough,
                          double  z_therm,
                          double  f_surf_layer,
                          double *KM_d,
                          double *KH_d,
                          int     num);
