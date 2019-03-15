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
// physics module interface
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

#include <string>

#include "config_file.h"
#include "define.h"
#include "esp.h"
#include "grid.h"
#include "log_writer.h"
#include "simulation_setup.h"
#include "storage.h"

// return name of module for storage to output files
std::string phy_modules_get_name();

// can print out configurations
void phy_modules_print_config();

// add config variables for this module to main config reader
bool phy_modules_generate_config(config_file& config_reader);

// allocate and free memory locally and on device for this module
// register arrays that need to be updated in the dynamical core RK update steps
bool phy_modules_init_mem(const ESP& esp, device_RK_array_manager& phy_modules_core_arrays);
bool phy_modules_free_mem();

// initialise data
bool phy_modules_init_data(const ESP& esp, const SimulationSetup& sim, storage* s);


// Dynamical core loop functions
bool phy_modules_dyn_core_loop_init(const ESP& esp);
bool phy_modules_dyn_core_loop_slow_modes(const ESP&             esp,
                                          const SimulationSetup& sim,
                                          int                    nstep, // Step number
                                          double                 times);                // Time-step [s]
bool phy_modules_dyn_core_loop_fast_modes(const ESP&             esp,
                                          const SimulationSetup& sim,
                                          int                    nstep, // Step number
                                          double                 time_step);            // Time-step [s]

bool phy_modules_dyn_core_loop_end(const ESP& esp);

// Physics loop
bool phy_modules_phy_loop(ESP&                   esp,
                          const SimulationSetup& sim,
                          int                    nstep, // Step number
                          double                 time_step);            // Time-step [s]


bool phy_modules_store(const ESP& esp, storage& s);

bool phy_modules_store_init(storage& s);
