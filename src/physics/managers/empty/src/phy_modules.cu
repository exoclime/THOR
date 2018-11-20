// empty physical modules template
// copy to your module folder and fill in the blanks.
// see examples in modules_simple_template
// and modules_complex_template

#include "phy_modules.h"

std::string phy_modules_get_name() {
    return std::string("empty");
}

void phy_modules_print_config() {
    printf("  empty physics module\n");
}

bool phy_modules_init_mem(const ESP&               esp,
                          device_RK_array_manager& phy_modules_core_arrays) {

    // initialise memory on device and host

    // register with
    // -> phy_modules_core_arrays(double* array, double* arrayk, double* arrayi, int     dimensions)
    // the device arrays that need to be updated by Runge Kutta scheme in main dynamical core loop
    return true;
}

bool phy_modules_init_data(const ESP&     esp,
                           const XPlanet& planet,
                           storage*       s) {

    // initialise the initial conditions

    // if storage * s pointer is not a nullptr, reload your initial conditions from that file

    return true;
}

bool phy_modules_generate_config(config_file& config_reader) {

    // register configuration variables with config_reader.
    // define variable to be written to, as global or as variale in a subclass
    // config_reader.append_config_var("name_in_config_file", <variable to store result to>,
    // <default value>);

    return true;
}

// ***********************************************************************************************
// Dynamical core update steps

bool phy_modules_dyn_core_loop_init(const ESP& esp) {

    // called before starting dynamical core update to swap buffers, zero out arrays, etc...

    return true;
}

// Slow modes, called after hyperdiffusion step
bool phy_modules_dyn_core_loop_slow_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep, // Step number
                                          double         times, // Time-step [s]
                                          bool           HyDiff) {

    return true;
}

// fast mode
bool phy_modules_dyn_core_loop_fast_modes(const ESP&     esp,
                                          const XPlanet& planet,
                                          int            nstep,     // Step number
                                          double         time_step) { // Time-step [s]

    return true;
}


bool phy_modules_dyn_core_loop_end(const ESP& esp) {

    // called after dynamical core to update buffers

    return true;
}

bool phy_modules_phy_loop(ESP&           esp,
                          const XPlanet& planet,
                          int            nstep,
                          double         time_step) {


    return true;
}

bool phy_modules_store_init(storage& s) {
    // store initial conditions to initial file
    return true;
}


bool phy_modules_store(const ESP& esp, storage& s) {
    // store state to output file
    return true;
}


bool phy_modules_free_mem() {
    // free memory at end of run
    return true;
}
