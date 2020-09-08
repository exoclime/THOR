// **********************************************************************************
//
// Example of external module to reuse phy module code at various places
// This pushes the module code in another file, with a standard structure, that make it easy to
// put modules in a list and reuse them

#include "phy_modules.h"

#include "boundary_layer.h"
#include "chemistry.h"
#include "radiative_transfer.h"

#include "log_writer.h"

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

bool       boundary_layer_enabled         = false;
const bool boundary_layer_enabled_default = false;

boundary_layer bl;

#ifdef HAS_ALFRODULL
#    include "two_streams_radiative_transfer.h"

bool       alfrodull_enabled         = false;
const bool alfrodull_enabled_default = false;

two_streams_radiative_transfer tsrt;
#endif // HAS_ALFRODULL

std::string phy_modules_get_name() {
    return std::string("multi");
}

void phy_modules_print_config() {
    log::printf("  multi physics module, with radiative transfer, chemistry, boundary layer\n");
    log::printf("   Radiative Transfer module: %s.\n",
                radiative_transfer_enabled ? "true" : "false");
    log::printf("   Boundary Layer module: %s.\n", boundary_layer_enabled ? "true" : "false");
    log::printf("\n");

    if (radiative_transfer_enabled)
        rt.print_config();

    log::printf("   Chemistry module: %s.\n", chemistry_enabled ? "true" : "false");
    if (chemistry_enabled)
        chem.print_config();
    if (boundary_layer_enabled)
        bl.print_config();

#ifdef HAS_ALFRODULL
    log::printf("   Two Stream Radiative Transfer module: %s.\n",
                alfrodull_enabled ? "true" : "false");
    if (alfrodull_enabled)
        tsrt.print_config();
#endif // HAS_ALFRODULL
}


bool phy_modules_init_mem(const ESP& esp, device_RK_array_manager& phy_modules_core_arrays) {
    // initialise all the modules memory

    bool out = true;

    if (radiative_transfer_enabled)
        rt.initialise_memory(esp, phy_modules_core_arrays);

    if (chemistry_enabled)
        chem.initialise_memory(esp, phy_modules_core_arrays);
    if (boundary_layer_enabled)
        bl.initialise_memory(esp, phy_modules_core_arrays);

#ifdef HAS_ALFRODULL
    if (alfrodull_enabled)
        tsrt.initialise_memory(esp, phy_modules_core_arrays);
#endif // HAS_ALFRODULL

    return out;
}

bool phy_modules_init_data(const ESP& esp, const SimulationSetup& sim, storage* s) {
    bool out = true;
    // initialise all the modules data

    // if (s != nullptr) {
    //     // load initialisation data from storage s
    // }

    if (radiative_transfer_enabled)
        out &= rt.initial_conditions(esp, sim, s);

    if (boundary_layer_enabled)
        out &= bl.initial_conditions(esp, sim, s);

    if (chemistry_enabled)
        out &= chem.initial_conditions(esp, sim, s);

#ifdef HAS_ALFRODULL
    if (alfrodull_enabled)
        out &= tsrt.initial_conditions(esp, sim, s);
#endif // HAS_ALFRODULL

    return out;
}

bool phy_modules_generate_config(config_file& config_reader) {
    bool out = true;

    config_reader.append_config_var(
        "radiative_transfer", radiative_transfer_enabled, radiative_transfer_enabled_default);

    rt.configure(config_reader);

    config_reader.append_config_var("chemistry", chemistry_enabled, chemistry_enabled_default);

    chem.configure(config_reader);

    config_reader.append_config_var(
        "boundary_layer", boundary_layer_enabled, boundary_layer_enabled_default);

    bl.configure(config_reader);


#ifdef HAS_ALFRODULL
    config_reader.append_config_var(
        "two_streams_radiative_transfer", alfrodull_enabled, alfrodull_enabled_default);

    tsrt.configure(config_reader);
#endif // HAS_ALFRODULL

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
                          kernel_diagnostics&    diag,
                          int                    nstep,
                          double                 time_step) {
    // run all the modules main loop
    bool out = true;

    if (chemistry_enabled)
        chem.phy_loop(esp, sim, diag, nstep, time_step);

    if (radiative_transfer_enabled)
        rt.phy_loop(esp, sim, diag, nstep, time_step);

#ifdef HAS_ALFRODULL
    if (alfrodull_enabled)
        tsrt.phy_loop(esp, sim, diag, nstep, time_step);
#endif // HAS_ALFRODULL

    if (boundary_layer_enabled)
        bl.phy_loop(esp, sim, diag, nstep, time_step);


    return out;
}

bool phy_modules_store_init(storage& s) {
    // radiative transfer option
    s.append_value(radiative_transfer_enabled ? 1.0 : 0.0,
                   "/radiative_transfer",
                   "-",
                   "Using radiative transfer");

    rt.store_init(s);

    // chemistry option
    s.append_value(chemistry_enabled ? 1.0 : 0.0, "/chemistry", "-", "Using relaxation chemistry");

    // store state of chemistry variables, so that we know it was disabled
    chem.store_init(s);

    s.append_value(
        boundary_layer_enabled ? 1.0 : 0.0, "/boundary_layer", "-", "Using boundary layer");

    bl.store_init(s);

#ifdef HAS_ALFRODULL

    s.append_value(alfrodull_enabled ? 1.0 : 0.0,
                   "/two_streams_radiative_transfer",
                   "-",
                   "Using two stream radiative transfer");

    tsrt.store_init(s);
#endif // HAS_ALFRODULL

    return true;
}

bool phy_modules_store(const ESP& esp, storage& s) {

    if (radiative_transfer_enabled)
        rt.store(esp, s);

    if (chemistry_enabled)
        chem.store(esp, s);

    if (boundary_layer_enabled)
        bl.store(esp, s);

#ifdef HAS_ALFRODULL
    if (alfrodull_enabled)
        tsrt.store(esp, s);
#endif // HAS_ALFRODULL

    return true;
}


bool phy_modules_free_mem() {
    // generate all the modules config
    bool out = true;

    if (radiative_transfer_enabled)
        rt.free_memory();

    if (chemistry_enabled)
        chem.free_memory();
    if (boundary_layer_enabled)
        bl.free_memory();

#ifdef HAS_ALFRODULL
    if (alfrodull_enabled)
        tsrt.free_memory();
#endif // HAS_ALFRODULL

    return out;
}
