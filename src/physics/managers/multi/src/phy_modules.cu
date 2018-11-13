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

std::string phy_modules_get_name()
{
    return std::string("multi");
}

void phy_modules_print_config()
{
    printf("  multi physics module, with radiative transfer\n");
    rt.print_config();
}


        
bool phy_modules_init_mem(const ESP& esp) {
    // initialise all the modules memory

    bool out = true;

    rt.initialise_memory(esp);

    return out;
}

bool phy_modules_init_data(const ESP&     esp,
                           const XPlanet& planet) {
    bool out = true;
    // initialise all the modules data

    rt.initial_conditions(esp, planet);

    return out;
}

bool phy_modules_generate_config(config_file& config_reader) {
    bool out = true;

    rt.configure(config_reader);

    return out;
}

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
) {
    // run all the modules main loop
    bool out = true;

    rt.loop(esp,
            nstep,
            core_benchmark,
            time_step,
            Omega,
            Cp,
            Rd,
            mu,
            kb,
            P_Ref,
            Gravit,
            A);

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
