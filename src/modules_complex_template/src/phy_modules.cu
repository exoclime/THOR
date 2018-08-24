// **********************************************************************************
// 
// Example of external module if we want to reuse phy module code at various places
// This pushes the module code in another file, with a standard structure, that make it easy to
// put modules in a list and reuse them

#include "phy_modules.h"

#include "profx_RT.h"

#include "radiative_transfer.h"

#include <math.h>
#include <vector>
#include <memory>

// define all the modules we want to use
radiative_transfer rt;

using std::vector;

vector<std::unique_ptr<phy_module_base>> modules;





bool phy_modules_init_mem(const ESP & esp)
{
    // initialise all the modules memory
    
    bool out = true;
    
    for (auto & m : modules)
    {
        out |= m->initialise_memory(esp);
    }
    
    return out;
}

bool phy_modules_init_data()
{
    // initialise all the modules data

    
    bool out = true;
    
    for (auto & m : modules)
    {
        out |= m->initial_conditions();
    }
    
    return out;
    
}

bool phy_modules_generate_config(config_file & config_reader)
{
    // This is our first entry point to modules, initialise stuff here
    modules.push_back( std::unique_ptr<phy_module_base>(new radiative_transfer()));
    
    // generate all the modules config
    bool out = true;
    
    for (auto & m : modules)
    {
        out |= m->configure(config_reader);
    }
    
    return out;
}

bool phy_modules_mainloop(ESP & esp,
                          int    nstep       , // Step number
                          int    hstest      , // Held-Suarez test option
                          double time_step   , // Time-step [s]
                          double Omega       , // Rotation rate [1/s]
                          double Cp          , // Specific heat capacity [J/kg/K]
                          double Rd          , // Gas constant [J/kg/K]
                          double Mmol        , // Mean molecular mass of dry air [kg]
                          double mu          , // Atomic mass unit [kg]
                          double kb          , // Boltzmann constant [J/K]
                          double P_Ref       , // Reference pressure [Pa]
                          double Gravit      , // Gravity [m/s^2]
                          double A           // Planet radius [m]
    )
{
    // run all the modules main loop
    bool out = true;
    
    for (auto & m : modules)
    {
        out |= m->loop(esp,
                       nstep       , 
                       hstest      , 
                       time_step   ,
                       Omega       ,
                       Cp          ,
                       Rd          ,
                       Mmol        ,
                       mu          ,
                       kb          ,
                       P_Ref       ,
                       Gravit      ,
                       A           );
    }
    
    return true;
}

bool phy_modules_store()
{

    return true;
}


bool phy_modules_free_mem()
{
    // generate all the modules config
    bool out = true;
    
    for (auto & m : modules)
    {
        out |= m->free_memory();
    }

    // clear all modules
    modules.clear();
    
    return out;
}
