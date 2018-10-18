// empty physical modules template
// copy to your module folder and fill in the blanks.
// see examples in modules_simple_template
// and modules_complex_template

#include "phy_modules.h"





bool phy_modules_init_mem(const ESP & esp)
{

    return true;
}

bool phy_modules_init_data(const ESP & esp,
                           const XPlanet & planet)
{

    return true;
}

bool phy_modules_generate_config(config_file & config_reader)
{

    return true;
}

bool phy_modules_mainloop(ESP & esp,
                          int    nstep       , // Step number
                          int core_benchmark      , // Held-Suarez test option
                          double time_step   , // Time-step [s]
                          double Omega       , // Rotation rate [1/s]
                          double Cp          , // Specific heat capacity [J/kg/K]
                          double Rd          , // Gas constant [J/kg/K]
                          double mu          , // Atomic mass unit [kg]
                          double kb          , // Boltzmann constant [J/K]
                          double P_Ref       , // Reference pressure [Pa]
                          double Gravit      , // Gravity [m/s^2]
                          double A           // Planet radius [m]
    )
{

    return true;
}

bool phy_modules_store_init(storage & s)
{
    return true;
}


bool phy_modules_store(const ESP & esp, storage & s)
{

    return true;
}


bool phy_modules_free_mem()
{
    return true;
}
