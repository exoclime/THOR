#pragma once

#include "esp.h"
#include "config_file.h"


class phy_module_base
{
public:
    phy_module_base() {};

    ~phy_module_base() {};


    virtual bool initialise_memory(const ESP & esp) = 0;
    virtual bool initial_conditions() = 0;

    // TBD, how does it get data? friend of ESP ? grid ?
    virtual bool loop(ESP & esp,
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
              double A           // Planet radius [m]);
        ) = 0; 

    // TBD
    virtual bool store() = 0; //? should be "get_data()"? "store_data(h5file)"?
    
    virtual bool  configure(config_file & config_reader) = 0;

    virtual bool free_memory() = 0;
};


    
