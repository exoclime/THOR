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
                  int                    nstep, // Step number
                  double                 time_step);            // Time-step [s]

    bool store(const ESP &esp, storage &s);

    bool store_init(storage &s);

    bool configure(config_file &config_reader);

    virtual bool free_memory();

    void print_config();

private:
    // Config options
    double Tstar_config            = 4520;  // Star effective temperature [K]
    double planet_star_dist_config = 0.015; // Planet-star distance [au]
    double radius_star_config      = 0.667; // Star radius [Rsun]
    double diff_ang_config         = 0.5;   // Diffusivity angle (1 / diffusivity factor): 0.5-1.0
    double Tlow_config =
        970; // Lower boundary temperature: upward flux coming from the planet's interior
    double albedo_config     = 0.18;   // Bond albedo
    double tausw_config      = 532.0;  // Absorption coefficient for the shortwaves
    double taulw_config      = 1064.0; // Absorption coefficient for the longwaves
    bool   latf_lw_config    = false;  // use sin^2(lat) dependence for lw opacity
    double taulw_pole_config = 1064.0; // Absorption coefficient for the longwave (poles)
    double n_lw_config       = 2.0;    // power law dependence for unmixed absorbers in LW
    double n_sw_config       = 1.0;    // power law dependence for mixed/unmixed absorbers in SW
    double f_lw_config       = 0.5;    // fraction of taulw in well-mixed absorber


    bool   sync_rot_config    = true;     // is planet syncronously rotating?
    double mean_motion_config = 1.991e-7; // orbital mean motion (rad/s)
    double true_long_i_config = 0;        // initial true longitude of planet (rad)
    double ecc_config         = 0;        // orbital eccentricity
    double obliquity_config   = 0;        // obliquity (tilt of spin axis) (rad)
    double alpha_i_config     = 0;        // initial right asc of host star (relative to long = 0)
    double longp_config       = 0;        // longitude of periastron (rad)

    bool   surface_config  = false; // use solid/liquid surface at altitude 0
    double Csurf_config    = 1e7;   // heat capacity of surface (J K^-1 m^-2)
    bool   rt1Dmode_config = false; // 1D mode=all columns are irradiated identically

    // Rad trans
    double Tstar            = 4520;  // Star effective temperature [K]
    double planet_star_dist = 0.015; // Planet-star distance [au]
    double radius_star      = 0.667; // Star radius [Rsun]
    double diff_ang         = 0.5;   // Diffusivity factor: 0.5-1.0
    double Tlow = 970; // Lower boundary temperature: upward flux coming from the planet's interior
    double albedo     = 0.18;   // Bond albedo
    double tausw      = 532.0;  // Absorption coefficient for the shortwaves
    double taulw      = 1064.0; // Absorption coefficient for the longwaves
    double taulw_pole = 1064.0;
    bool   latf_lw    = false;
    double n_lw       = 2.0; // power law dependence for unmixed absorbers in LW
    double n_sw       = 1.0; // power law dependence for mixed/unmixed absorbers in SW
    double f_lw       = 0.5; // fraction of taulw in well-mixed absorber

    bool    rt1Dmode;
    bool    surface;
    double  Csurf;
    double *surf_flux_d;
    double *Tsurface_d;
    double *Tsurface_h;

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

    // orbit/insolation properties
    bool   sync_rot       = true;     // is planet syncronously rotating?
    double mean_motion    = 1.991e-7; // orbital mean motion (rad/s)
    double mean_anomaly_i = 0;        // initial mean anomaly at start (rad)
    double mean_anomaly   = 0;        // current mean anomaly of planet (rad)
    double true_long_i    = 0;        // initial true longitude of planet (rad)
    double ecc            = 0;        // orbital eccentricity
    double obliquity      = 0;        // obliquity (tilt of spin axis) (rad)
    double r_orb          = 1;        // orbital distance/semi-major axis
    double sin_decl       = 0;        // declination of host star (relative to equator)
    double cos_decl       = 1;
    double alpha_i        = 0; // initial right asc of host star (relative to long = 0)
    double alpha          = 0; // right asc of host star (relative to long = 0)
    double longp          = 0; // longitude of periastron (rad)

    double *insol_h;
    double *insol_d;
    double *insol_ann_h;
    double *insol_ann_d;

    //  These arrays are for temporary usage in RT code
    double *dtemp;
    double *phtemp;
    double *ttemp;
    double *thtemp;
    void    RTSetup(double Tstar_,
                    double planet_star_dist_,
                    double radius_star_,
                    double diff_ang_,
                    double Tlow_,
                    double albedo_,
                    double tausw_,
                    double taulw_,
                    bool   latf_lw_,
                    double taulw_pole_,
                    double n_lw_,
                    double n_sw_,
                    double f_lw_,
                    bool   sync_rot_,
                    double mean_motion_,
                    double true_long_i_,
                    double longp_,
                    double ecc_,
                    double alpha_i_,
                    double obliquity_,
                    double Omega,
                    bool   surface,
                    double Csurf,
                    bool   rt1Dmode,
                    double Tmean,
                    int    point_num);

    void update_spin_orbit(double time, double Omega);
};
