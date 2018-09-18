#pragma once

#include "phy_module_base.h"

class radiative_transfer : public phy_module_base
{
public:
    radiative_transfer();
    ~radiative_transfer();

    bool initialise_memory(const ESP & esp);
    bool initial_conditions();

    // TBD, how does it get data? friend of ESP ? grid ?
    bool loop(ESP & esp,
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
        );

    bool store(storage & s);

    bool store_init(storage & s);

    bool configure(config_file & config_reader);

    virtual bool free_memory();

private:


    // Rad Trans options
    double Tstar            = 4520;       // Star effective temperature [K]
    double planet_star_dist = 0.015;      // Planet-star distance [au]
    double radius_star      = 0.667;      // Star radius [Rsun]
    double diff_fac         = 0.5;        // Diffusivity factor: 0.5-1.0
    double Tlow             = 970;        // Lower boundary temperature: upward flux coming from the planet's interior
    double albedo           = 0.18;       // Bond albedo
    double tausw            = 532.0;      // Absorption coefficient for the shortwaves
    double taulw            = 1064.0;     // Absorption coefficient for the longwaves

    // double resc_flx       ; TODO: check if this is still used
    double incflx         ;
    //  Arrays used in RT code
    double *fnet_up_d     ;
    double *fnet_dn_d     ;
    double *tau_d         ;

    //  These arrays are for temporary usage in RT code
    double *dtemp         ;
    double *phtemp        ;
    double *ttemp         ;
    double *thtemp        ;
    void RTSetup(double Tstar_           ,
                 double planet_star_dist_,
                 double radius_star_     ,
                 double diff_fac_        ,
                 double Tlow_            ,
                 double albedo_          ,
                 double tausw_           ,
                 double taulw_           );
};
