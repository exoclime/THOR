// **********************************************************************************
//
// Example of external module to reuse phy module code at various places
// This pushes the module code in another file, with a standard structure, that make it easy to
// put modules in a list and reuse them

#include "phy_modules.h"

#include "chemistry.h"
#include "radiative_transfer.h"
#include "log_writer.h"
#include "radiative_transfer.h"

#include <math.h>
#include <memory>
#include <vector>

// define all the modules we want to use

bool       radiative_transfer_enabled         = false;
const bool radiative_transfer_enabled_default = false;

radiative_transfer rt;

bool       chemistry_enabled         = false;
const bool chemistry_enabled_default = false;

chemistry chem;

std::string phy_modules_get_name() {
    return std::string("multi");
}

void phy_modules_print_config() {
    log::printf("  multi physics module, with radiative transfer and chemistry\n");
    log::printf("   Radiative Transfer module: %s.\n", radiative_transfer_enabled ? "true" : "false");

    if (radiative_transfer_enabled)
        rt.print_config();

    log::printf("   Chemistry module: %s.\n", chemistry_enabled ? "true" : "false");
    if (chemistry_enabled)
    	chem.print_config();
}


bool phy_modules_init_mem(const ESP&               esp,
                          device_RK_array_manager& phy_modules_core_arrays) {
    // initialise all the modules memory

    bool out = true;

    if (radiative_transfer_enabled)
        rt.initialise_memory(esp, phy_modules_core_arrays);

    if (chemistry_enabled)
        chem.initialise_memory(esp, phy_modules_core_arrays);


    return out;
}

bool phy_modules_init_data(const ESP&             esp,
                           const SimulationSetup& sim,
                           storage*               s) {
    bool out = true;
    // initialise all the modules data

    if (s != nullptr) {
        // load initialisation data from storage s
    }

    if (radiative_transfer_enabled)
        out &= rt.initial_conditions(esp, sim);

    if (chemistry_enabled)
        out &= chem.initial_conditions(esp, sim);
    return out;
}

bool phy_modules_generate_config(config_file& config_reader) {
    bool out = true;

    config_reader.append_config_var("radiative_transfer", radiative_transfer_enabled, radiative_transfer_enabled_default);

    rt.configure(config_reader);

    config_reader.append_config_var("chemistry", chemistry_enabled, chemistry_enabled_default);

    chem.configure(config_reader);

    return out;
}

bool phy_modules_dyn_core_loop_init(const ESP& esp) {

    if (chemistry_enabled)
        chem.dyn_core_loop_init(esp);

    return true;
}

bool phy_modules_dyn_core_loop_slow_modes(const ESP&             esp,
                                          const SimulationSetup& sim,
                                          int                    nstep, // Step number
                                          double                 times) {               // Time-step [s]

    if (chemistry_enabled)
        chem.dyn_core_loop_slow_modes(esp, sim, nstep, times);

    return true;
}

bool phy_modules_dyn_core_loop_fast_modes(const ESP&             esp,
                                          const SimulationSetup& sim,
                                          int                    nstep, // Step number
                                          double                 time_step) {           // Time-step [s]

    if (chemistry_enabled)
        chem.dyn_core_loop_fast_modes(esp, sim, nstep, time_step);

    return true;
}

bool phy_modules_dyn_core_loop_end(const ESP& esp) {

    if (chemistry_enabled)
        chem.dyn_core_loop_end(esp);

    return true;
}


bool phy_modules_phy_loop(ESP&                   esp,
                          const SimulationSetup& sim,
                          int                    nstep,
                          double                 time_step) {
    // run all the modules main loop
    bool out = true;

    if (chemistry_enabled)
        chem.phy_loop(esp, sim, nstep, time_step);

    if (radiative_transfer_enabled)
        rt.phy_loop(esp, sim, nstep, time_step);

    return out;
}

bool phy_modules_store_init(storage& s) {
    // radiative transfer option
    s.append_value(radiative_transfer_enabled ? 1.0 : 0.0, "/radiativetransfer", "-", "Using radiative transfer");

    rt.store_init(s);

    // chemistry option
    s.append_value(chemistry_enabled ? 1.0 : 0.0, "/chemistry", "-", "Using relaxation chemistry");

    // store state of chemistry variables, so that we know it was disabled
    chem.store_init(s);

    return true;
}

bool phy_modules_store(const ESP& esp, storage& s) {

    if (radiative_transfer_enabled)
        rt.store(esp, s);

    if (chemistry_enabled)
        chem.store(esp, s);

    return true;
}


bool phy_modules_free_mem() {
    // generate all the modules config
    bool out = true;

    if (radiative_transfer_enabled)
        rt.free_memory();

    if (chemistry_enabled)
        chem.free_memory();


    return out;
}
