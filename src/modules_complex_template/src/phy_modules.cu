#include "phy_modules.h"
#include <math.h>

// Rad Trans options
double Tstar            = 4520;       // Star effective temperature [K]
double planet_star_dist = 0.015;      // Planet-star distance [au]
double radius_star      = 0.667;      // Star radius [Rsun]
double diff_fac         = 0.5;        // Diffusivity factor: 0.5-1.0
double Tlow             = 970;        // Lower boundary temperature: upward flux coming from the planet's interior
double albedo           = 0.18;       // Bond albedo
double tausw            = 532.0;      // Absorption coefficient for the shortwaves
double taulw            = 1064.0;     // Absorption coefficient for the longwaves

//  Arrays used in RT code
double *fnet_up_d     ;
double *fnet_dn_d     ;
double *tau_d         ;

//  These arrays are for temporary usage in RT code
double *dtemp         ;
double *phtemp        ;
double *ttemp         ;
double *thtemp        ;

RTSetup(double Tstar_           ,
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


bool phy_modules_init_mem()
{
//  Rad Transfer
    cudaMalloc((void **)&fnet_up_d   , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&fnet_dn_d   , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&tau_d       , nv * point_num * 2 *  sizeof(double));

    cudaMalloc((void **)&phtemp      , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&thtemp      , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&ttemp       , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&dtemp       , nv * point_num *     sizeof(double));

    
    return true;
}

bool phy_module_init_data()
{

    return true;
}

bool phy_modules_generate_config(config_file & config_reader)
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

bool phy_modules_mainloop()
{

    return true;
}

bool phy_modules_store()
{

    return true;
}


bool phy_modules_free_mem()
{
    cudaFree(fnet_up_d);
    cudaFree(fnet_dn_d);
    cudaFree(tau_d);

    cudaFree(phtemp);
    cudaFree(thtemp);
    cudaFree(ttemp );
    cudaFree(dtemp);

    
}
