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

// allocate and free memory localy and on device for this module
bool phy_modules_init_mem(const ESP& esp);
bool phy_modules_free_mem();

// initialise data
bool phy_modules_init_data(const ESP& esp,
                           const XPlanet& planet,
                           storage * s = nullptr);


// TODO: what do we need here? esp,planet, grid???
bool phy_modules_mainloop(ESP&            esp,
                          int             nstep,          // Step number
                          benchmark_types core_benchmark, // Held-Suarez test option
                          double          time_step,      // Time-step [s]
                          double          Omega,          // Rotation rate [1/s]
                          double          Cp,             // Specific heat capacity [J/kg/K]
                          double          Rd,             // Gas constant [J/kg/K]
                          double          mu,             // Atomic mass unit [kg]
                          double          kb,             // Boltzmann constant [J/K]
                          double          P_Ref,          // Reference pressure [Pa]
                          double          Gravit,         // Gravity [m/s^2]
                          double          A               // Planet radius [m]
);


bool phy_modules_store(const ESP& esp,
                       storage&   s);

bool phy_modules_store_init(storage& s);
