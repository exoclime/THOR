#pragma once

#include "config_file.h"
#include "define.h"
#include "dyn/phy_modules_device.h"
#include "esp.h"
#include "log_writer.h"
#include "simulation_setup.h"
#include "storage.h"

class phy_module_base
{
public:
    phy_module_base(){};

    ~phy_module_base(){};


    virtual bool initialise_memory(const ESP&               esp,
                                   device_RK_array_manager& phy_modules_core_arrays)        = 0;
    virtual bool initial_conditions(const ESP& esp, const SimulationSetup& sim, storage* s) = 0;

    virtual bool dyn_core_loop_init(const ESP& esp) {
        return true;
    };
    virtual bool dyn_core_loop_slow_modes(const ESP&             esp,
                                          const SimulationSetup& sim,

                                          int    nstep, // Step number
                                          double times) // Time-step [s]
    {
        return true;
    };
    virtual bool dyn_core_loop_fast_modes(const ESP&             esp,
                                          const SimulationSetup& sim,
                                          int                    nstep, // Step number
                                          double                 time_step)             // Time-step [s]
    {
        return true;
    };
    virtual bool dyn_core_loop_end(const ESP& esp) {
        return true;
    };

    virtual bool phy_loop(ESP&                   esp,
                          const SimulationSetup& sim,
                          int                    nstep,    // Step number
                          double                 time_step // Time-step [s]
                          ) = 0;


    virtual bool store(const ESP& esp, storage& s) = 0;

    virtual bool configure(config_file& config_reader) = 0;

    virtual bool free_memory() = 0;
};
