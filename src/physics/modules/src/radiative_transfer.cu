#include "radiative_transfer.h"

#include "profx_RT.h"

radiative_transfer::radiative_transfer()
{


}

radiative_transfer::~radiative_transfer()
{


}

bool radiative_transfer::initialise_memory(const ESP & esp)
{
    //  Rad Transfer
    cudaMalloc((void **)&fnet_up_d   , esp.nvi * esp.point_num *     sizeof(double));
    cudaMalloc((void **)&fnet_dn_d   , esp.nvi * esp.point_num *     sizeof(double));
    cudaMalloc((void **)&tau_d       , esp.nv * esp.point_num * 2 *  sizeof(double));

    cudaMalloc((void **)&phtemp      , esp.nvi * esp.point_num *     sizeof(double));
    cudaMalloc((void **)&thtemp      , esp.nvi * esp.point_num *     sizeof(double));
    cudaMalloc((void **)&ttemp       , esp.nv * esp.point_num *     sizeof(double));
    cudaMalloc((void **)&dtemp       , esp.nv * esp.point_num *     sizeof(double));

    
    return true;
}


bool radiative_transfer::free_memory()
{
    cudaFree(fnet_up_d);
    cudaFree(fnet_dn_d);
    cudaFree(tau_d);
    
    cudaFree(phtemp);
    cudaFree(thtemp);
    cudaFree(ttemp );
    cudaFree(dtemp);

    return true;
}

bool radiative_transfer::initial_conditions()
{
    RTSetup(Tstar            ,
            planet_star_dist ,
            radius_star      ,
            diff_fac         ,
            Tlow             ,
            albedo           ,
            tausw            ,
            taulw            );
    return true;

}

bool radiative_transfer::loop(ESP & esp,
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
    )
{
//
//  Number of threads per block.
    const int NTH = 256;

//  Specify the block sizes.
    dim3 NB((esp.point_num / NTH) + 1, esp.nv, 1);
    dim3 NBRT((esp.point_num/NTH) + 1, 1, 1);
    
    rtm_dual_band <<< NBRT, NTH >>> (esp.pressure_d   ,
      //rtm_dual_band <<< 1,1 >>> (pressure_d         ,
                                     esp.Rho_d        ,
                                     esp.temperature_d,
                                     fnet_up_d    ,
                                     fnet_dn_d    ,
                                     tau_d        ,
                                     Gravit       ,
                                     Cp           ,
                                     esp.lonlat_d     ,
                                     esp.Altitude_d   ,
                                     esp.Altitudeh_d  ,
                                     phtemp       ,
                                     dtemp        ,
                                     ttemp        ,
                                     thtemp       ,
                                     time_step    ,
                                     Tstar        ,
                                     planet_star_dist,
                                     radius_star  ,
                                     diff_fac     ,
                                     Tlow         ,
                                     albedo       ,
                                     tausw        ,
                                     taulw        ,
                                     incflx       ,
                                     P_Ref        ,
                                     esp.point_num    ,
                                     esp.nv           ,
                                     esp.nvi          ,
                                     A             );


    return true;
}

bool radiative_transfer::configure(config_file & config_reader)
{
config_reader.append_config_var("Tstar", Tstar, Tstar);
    config_reader.append_config_var("planet_star_dist", planet_star_dist, planet_star_dist);
    config_reader.append_config_var("radius_star", radius_star, radius_star);
    config_reader.append_config_var("diff_fac", diff_fac, diff_fac);
    config_reader.append_config_var("Tlow", Tlow, Tlow);
    config_reader.append_config_var("albedo", albedo, albedo);
    config_reader.append_config_var("tausw", tausw, tausw);
    config_reader.append_config_var("taulw", taulw, taulw);

    return true;
}

bool radiative_transfer::store(storage & s)
{
    
    return true;
}

void radiative_transfer::RTSetup(double Tstar_           ,
                                 double planet_star_dist_,
                                 double radius_star_     ,
                                 double diff_fac_        ,
                                 double Tlow_            ,
                                 double albedo_          ,
                                 double tausw_           ,
                                 double taulw_           ) {
    
    double bc = 5.677036E-8; // Stefan–Boltzmann constant [W m−2 K−4]
    
    Tstar = Tstar_;
    planet_star_dist = planet_star_dist_*149597870.7;
    radius_star = radius_star_*695508;
    diff_fac = diff_fac_;
    Tlow = Tlow_;
    albedo = albedo_;
    tausw = tausw_;
    taulw = taulw_;
    double resc_flx = pow(radius_star/planet_star_dist,2.0);
    incflx = resc_flx*bc*Tstar*Tstar*Tstar*Tstar;
}
