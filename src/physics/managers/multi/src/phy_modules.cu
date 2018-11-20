// **********************************************************************************
//
// Example of external module to reuse phy module code at various places
// This pushes the module code in another file, with a standard structure, that make it easy to
// put modules in a list and reuse them

#include "phy_modules.h"

#include "radiative_transfer.h"

#include <math.h>
#include <memory>
#include <vector>

// define all the modules we want to use
radiative_transfer rt;

std::string phy_modules_get_name() {
    return std::string("multi");
}

void phy_modules_print_config() {
    printf("  multi physics module, with radiative transfer and chemistry\n");
    rt.print_config();
}


bool phy_modules_init_mem(const ESP&               esp,
                          device_RK_array_manager& phy_modules_core_arrays) {
    // initialise all the modules memory

    bool out = true;

    rt.initialise_memory(esp, phy_modules_core_arrays);

    return out;
}

bool phy_modules_init_data(const ESP&     esp,
                           const XPlanet& planet,
                           storage*       s) {
    bool out = true;
    // initialise all the modules data

    if (s != nullptr) {
        // load initialisation data from storage s
    }


    out &= rt.initial_conditions(esp, planet);

    return out;
}

bool phy_modules_generate_config(config_file& config_reader) {
    bool out = true;

    rt.configure(config_reader);

    return out;
}

bool phy_modules_dyn_core_loop_init(const ESP& esp) {

    return true;
}

bool phy_modules_dyn_core_loop_slow_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep, // Step number
                                          double         times, // Time-step [s]
                                          double         mu,    // Atomic mass unit [kg]
                                          double         kb,    // Boltzmann constant [J/K]
                                          bool           HyDiff) {

    return true;
}

bool phy_modules_dyn_core_loop_fast_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep,     // Step number
                                          double         time_step, // Time-step [s]
                                          double         mu,        // Atomic mass unit [kg]
                                          double         kb) {

    return true;
}

bool phy_modules_dyn_core_loop_end(const ESP& esp) {

    return true;
}


bool phy_modules_phy_loop(ESP&           esp,
                          const XPlanet& planet,
                          int            nstep,
                          double         time_step,
                          double         mu,
                          double         kb


) {
    // run all the modules main loop
    bool out = true;
    
    rt.phy_loop(esp, planet, nstep, time_step, mu, kb);

    return out;
}

bool phy_modules_store_init(storage& s) {
    rt.store_init(s);

    return true;
}

bool phy_modules_store(const ESP& esp, storage& s) {

    rt.store(esp, s);

    return true;
}


bool phy_modules_free_mem() {
    // generate all the modules config
    bool out = true;

    rt.free_memory();

    return out;
}
