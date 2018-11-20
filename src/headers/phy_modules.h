// physical module interface

#pragma once

#include <string>

#include "config_file.h"
#include "define.h"
#include "esp.h"
#include "grid.h"
#include "planet.h"
#include "storage.h"

// return name of module for storage to output files
std::string phy_modules_get_name();

// can print out configurations
void phy_modules_print_config();

// add config variables for this module to main config reader
bool phy_modules_generate_config(config_file& config_reader);

// allocate and free memory locally and on device for this module
// register arrays that need to be updated in the dynamical core RK update steps
bool phy_modules_init_mem(const ESP&               esp,
                          device_RK_array_manager& phy_modules_core_arrays);
bool phy_modules_free_mem();

// initialise data
bool phy_modules_init_data(const ESP&     esp,
                           const XPlanet& planet,
                           storage*       s);


// Dynamical core loop functions
bool phy_modules_dyn_core_loop_init(const ESP& esp);
bool phy_modules_dyn_core_loop_slow_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep, // Step number
                                          double         times, // Time-step [s]
                                          bool           HyDiff);
bool phy_modules_dyn_core_loop_fast_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep, // Step number
                                          double         time_step);    // Time-step [s]

bool phy_modules_dyn_core_loop_end(const ESP& esp);

// Physics loop
// TODO: what do we need here? esp,planet, grid???

bool phy_modules_phy_loop(ESP&           esp,
                          const XPlanet& planet,
                          int            nstep, // Step number
                          double         time_step);    // Time-step [s]


bool phy_modules_store(const ESP& esp,
                       storage&   s);

bool phy_modules_store_init(storage& s);
