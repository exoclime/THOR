#pragma once

#include "phy_module_base.h"

class radiative_transfer : public phy_module_base
{
public:
    radiative_transfer();
    ~radiative_transfer();

    bool initialise_memory(const ESP & esp);
    bool initial_conditions(const ESP & esp,
                            const XPlanet & planet);

    // TBD, how does it get data? friend of ESP ? grid ?
    bool loop(ESP & esp,
              int    nstep       , // Step number
              int    hstest      , // Held-Suarez test option
              double time_step   , // Time-step [s]
              double Omega       , // Rotation rate [1/s]
              double Cp          , // Specific heat capacity [J/kg/K]
              double Rd          , // Gas constant [J/kg/K]
              double mu          , // Atomic mass unit [kg]
              double kb          , // Boltzmann constant [J/K]
              double P_Ref       , // Reference pressure [Pa]
              double Gravit      , // Gravity [m/s^2]
              double A           // Planet radius [m]);
        );

    bool store(const ESP & esp,
               storage & s);

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

    double *fnet_up_h     ;
    double *fnet_dn_h     ;
    double *tau_h         ;


    // orbit/insolation properties
    bool   sync_rot       = true     ;  // is planet syncronously rotating?
    double mean_motion    = 1.991e-7 ;  // orbital mean motion (rad/s)
    double mean_anomaly_i = 0        ;  // initial mean anomaly at start (rad)
    double mean_anomaly   = 0        ;  // current mean anomaly of planet (rad)
    double true_long_i    = 0        ;  // initial true longitude of planet (rad)
    double ecc            = 0        ;  // orbital eccentricity
    double obliquity      = 0        ;  // obliquity (tilt of spin axis) (rad)
    double r_orb          = 1        ;  // orbital distance/semi-major axis
    double sin_decl       = 0        ;  // declination of host star (relative to equator)
    double cos_decl       = 1        ;
    double alpha_i        = 0        ;  // initial right asc of host star (relative to long = 0)
    double alpha          = 0        ;  // right asc of host star (relative to long = 0)
    double longp          = 0        ;  // longitude of periastron (rad)

    double *insol_h;
    double *insol_d;

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
                 double taulw_           ,
                 bool   sync_rot_        ,
                 double mean_motion_     ,
                 double true_long_i_     ,
                 double longp_           ,
                 double ecc_             ,
                 double alpha_i_         ,
                 double obliquity_       ,
                 double Omega            ,
                 int    point_num);

    void update_spin_orbit(double time,
                           double Omega );
};
