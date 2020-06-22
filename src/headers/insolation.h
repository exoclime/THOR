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
// Description: Insolation computation
//
//
//
// Known limitations: None
//
// Known issues: None.
//
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
#pragma once

#include "cuda_device_memory.h"

#include "config_file.h"
#include "define.h"
#include "dyn/phy_modules_device.h"
#include "log_writer.h"
#include "simulation_setup.h"
#include "storage.h"

class ESP;

class Insolation
{
public:
    Insolation();

    // Called by any module to enable insolation computation
    void set_require() {
        enabled = true;
    };

    bool configure(config_file &config_reader);

    void print_config();

    bool initialise_memory(const ESP &esp, device_RK_array_manager &phy_modules_core_arrays);
    bool initial_conditions(const ESP &esp, const SimulationSetup &sim, storage *s);

    bool store_init(storage &s);

    bool phy_loop(ESP &                  esp,
                  const SimulationSetup &sim,
                  int                    nstep, // Step number
                  double                 time_step);            // Time-step [s]

    bool store(const ESP &esp, storage &s);

    // get pointer to data on device
    double *get_device_cos_zenith_angles() {
        return *cos_zenith_angles;
    };

    // fetch data to host and get pointer to data
    std::shared_ptr<double[]> get_host_cos_zenith_angles() {
        return cos_zenith_angles.get_host_data();
    }

    double get_r_orb() {
        return r_orb;
    }

    double get_mean_motion() {
        return mean_motion;
    }

private:
    bool enabled = false;


    // insolation computation vars from rt module
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

    bool   sync_rot_config    = true;     // is planet syncronously rotating?
    double mean_motion_config = 1.991e-7; // orbital mean motion (rad/s)
    double true_long_i_config = 0;        // initial true longitude of planet (rad)
    double ecc_config         = 0;        // orbital eccentricity
    double obliquity_config   = 0;        // obliquity (tilt of spin axis) (rad)
    double alpha_i_config     = 0;        // initial right asc of host star (relative to long = 0)
    double longp_config       = 0;        // longitude of periastron (rad)

    cuda_device_memory<double> cos_zenith_angles;

    void update_spin_orbit(double time, double Omega);
};
