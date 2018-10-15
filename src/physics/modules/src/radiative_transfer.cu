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

bool radiative_transfer::initial_conditions(const ESP & esp,
                                            const XPlanet & planet)
{
    RTSetup(Tstar            ,
            planet_star_dist ,
            radius_star      ,
            diff_fac         ,
            Tlow             ,
            albedo           ,
            tausw            ,
            taulw            ,
            sync_rot         ,
            mean_motion      ,
            true_long_i      ,
            longp            ,
            ecc              ,
            alpha_i          ,
            obliquity        ,
            planet.Omega     ,
            esp.point_num);
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

//  update global insolation properties if necessary
    if (sync_rot) {
      if (ecc > 1e-10) {
        update_spin_orbit(nstep*time_step, Omega);
      }
    } else {
      update_spin_orbit(nstep*time_step, Omega);
    }

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
                                     A             ,
                                     r_orb         ,
                                     alpha    ,  //current RA of star (relative to zero long on planet)
                                     alpha_i  ,
                                     sin_decl ,  //declination of star
                                     cos_decl ,
                                     sync_rot ,
                                     ecc      ,
                                     obliquity,
                                     insol_d);

    return true;
}

bool radiative_transfer::configure(config_file & config_reader)
{
    // basic star-planet properties
    config_reader.append_config_var("Tstar", Tstar, Tstar);
    config_reader.append_config_var("planet_star_dist", planet_star_dist, planet_star_dist);
    config_reader.append_config_var("radius_star", radius_star, radius_star);
    config_reader.append_config_var("diff_fac", diff_fac, diff_fac);
    config_reader.append_config_var("Tlow", Tlow, Tlow);
    config_reader.append_config_var("albedo", albedo, albedo);
    config_reader.append_config_var("tausw", tausw, tausw);
    config_reader.append_config_var("taulw", taulw, taulw);

    // orbit/insolation properties
    config_reader.append_config_var("sync_rot", sync_rot, sync_rot);
    config_reader.append_config_var("mean_motion", mean_motion, mean_motion);
    config_reader.append_config_var("alpha_i", alpha_i, alpha_i);
    config_reader.append_config_var("true_long_i", true_long_i, true_long_i);
    config_reader.append_config_var("ecc", ecc, ecc);
    config_reader.append_config_var("obliquity", obliquity, obliquity);
    config_reader.append_config_var("longp", longp, longp);

    return true;
}

bool radiative_transfer::store(const ESP & esp,
                               storage & s)
{
    cudaMemcpy(insol_h  , insol_d   , esp.point_num * sizeof(double), cudaMemcpyDeviceToHost);
    s.append_table( insol_h,
                     esp.point_num,
                     "/insol",
                     "W m^-2 s^-1",
                     "insolation (instantaneous)");

    return true;
}

bool radiative_transfer::store_init(storage & s)
{
    s.append_value(Tstar, "/Tstar", "K", "Temperature of host star");
    s.append_value(Tlow, "/Tlow", "K", "Temperature of interior heat flux");
    s.append_value(planet_star_dist/149597870.7, "/planet_star_dist", "au", "distance b/w host star and planet");
    s.append_value(radius_star/695508, "/radius_star", "R_sun", "radius of host star");
    s.append_value(diff_fac, "/diff_fac", "-", "diffusivity factor");
    s.append_value(albedo, "/albedo", "-", "bond albedo of planet");
    s.append_value(tausw, "/tausw", "-", "shortwave optical depth of deepest layer");
    s.append_value(taulw, "/taulw", "-", "longwave optical depth of deepest layer");
    s.append_value(sync_rot?1.0:0.0, "/sync_rot", "-", "enforce synchronous rotation");
    s.append_value(alpha_i*180/M_PI, "/alpha_i", "deg", "initial RA of host star");
    s.append_value(true_long_i*180/M_PI, "/true_long_i", "deg", "initial orbital position of planet");
    s.append_value(ecc, "/ecc", "-", "orbital eccentricity");
    s.append_value(obliquity*180/M_PI, "/obliquity", "deg", "tilt of spin axis");
    s.append_value(longp*180/M_PI, "/longp", "deg", "longitude of periastron");

    return true;
}

void radiative_transfer::RTSetup(double Tstar_           ,
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
                                 int    point_num        )
                               {

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

    sync_rot = sync_rot_;
    if (sync_rot) {
      mean_motion = Omega; //just set for the sync_rot, obl != 0 case
    } else {
      mean_motion = mean_motion_;
    }
    true_long_i = true_long_i_*M_PI/180.0;
    longp = longp_*M_PI/180.0;
    ecc = ecc_;
    double true_anomaly_i = fmod(true_long_i-longp,(2*M_PI));
    double ecc_anomaly_i = true2ecc_anomaly(true_anomaly_i,ecc);
    mean_anomaly_i = fmod(ecc_anomaly_i - ecc*sin(ecc_anomaly_i),(2*M_PI));
    alpha_i = alpha_i_*M_PI/180.0;
    obliquity = obliquity_*M_PI/180.0;

    insol_h        = (double*)malloc(point_num   * sizeof(double));
    cudaMalloc((void **)&insol_d, point_num *     sizeof(double));
}

void radiative_transfer::update_spin_orbit(double time  ,
                                           double Omega ) {

// Update the insolation related parameters for spin and orbit
  double ecc_anomaly, true_long;

  mean_anomaly = fmod((mean_motion*time + mean_anomaly_i),(2*M_PI));

  ecc_anomaly = fmod(solve_kepler(mean_anomaly, ecc),(2*M_PI));

  r_orb = calc_r_orb(ecc_anomaly, ecc);
  true_long = fmod((ecc2true_anomaly(ecc_anomaly, ecc) + longp),(2*M_PI));

  sin_decl = sin(obliquity)*sin(true_long);
  cos_decl = sqrt(1.0 - sin_decl*sin_decl);
  alpha = -Omega*time + true_long - true_long_i + alpha_i;
}
