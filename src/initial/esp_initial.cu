// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Build the class ESP (Exoclimes Simulation Platform)
//
//
// Description:
//   Declare and initialize variables in the model
//
// Method: -
//
//
// Known limitations: None.
//
//
// Known issues: None.
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "../headers/esp.h"
#include "hdf5.h"
#include <stdio.h>

__host__ ESP::ESP(int *point_local_    ,
                  int *maps_           ,
                  double *lonlat_      ,
                  double *Altitude_    ,
                  double *Altitudeh_   ,
                  double *nvecoa_      ,
                  double *nvecti_      ,
                  double *nvecte_      ,
                  double *areasT_      ,
                  double *areasTr_     ,
                  double *div_         ,
                  double *grad_        ,
                  double *func_r_      ,
                  int nl_region_       ,
                  int nr_              ,
                  int nv_              ,
                  int nvi_             ,
                  int nlat_            ,
                  int point_num_    ): nl_region(nl_region_), nr(nr_), point_num(point_num_), nv(nv_), nvi(nvi_), nlat(nlat_){

    point_local_h = point_local_;
    maps_h        = maps_       ;

    lonlat_h = lonlat_;

    Altitude_h = Altitude_ ;
    Altitudeh_h= Altitudeh_;

    nvecoa_h= nvecoa_    ;
    nvecti_h= nvecti_    ;
    nvecte_h= nvecte_    ;
    areasTr_h = areasTr_ ;
    areasT_h= areasT_    ;

    div_h = div_ ;
    grad_h= grad_;

    func_r_h = func_r_ ;

//
//  Allocate Data
    AllocData();

}

__host__ void ESP::AllocData(){


//
//  Description:
//
//  Allocate data on host and device.
//
//  Allocate data in host
//  Diagnostics
    Rho_h        = (double*)malloc(nv*point_num   * sizeof(double));
    pressure_h   = (double*)malloc(nv*point_num   * sizeof(double));
    temperature_h= (double*)malloc(nv*point_num   * sizeof(double));
    Mh_h         = (double*)malloc(nv*point_num*3 * sizeof(double));
    W_h          = (double*)malloc(nv*point_num   * sizeof(double));
    Wh_h         = (double*)malloc(nvi*point_num  * sizeof(double));

//  Allocate data in device
//  Grid
    cudaMalloc((void **)&point_local_d, 6 * point_num * sizeof(int));
    cudaMalloc((void **)&maps_d, (nl_region + 2)*(nl_region + 2)*nr * sizeof(int));

//  Operators
    cudaMalloc((void **)&nvecoa_d , 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecti_d , 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&nvecte_d , 6 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&areasTr_d, 6 * point_num * sizeof(double));
    cudaMalloc((void **)&func_r_d  , 3 * point_num * sizeof(double));
    cudaMalloc((void **)&div_d, 7 * 3 * point_num * sizeof(double));
    cudaMalloc((void **)&grad_d,7 * 3 * point_num * sizeof(double));

//  Altitude (grid)
    cudaMalloc((void **)&Altitude_d  , nv   * sizeof(double));
    cudaMalloc((void **)&Altitudeh_d , nvi  * sizeof(double));

//  Longitude-latitude
    cudaMalloc((void **)&lonlat_d  , 2 * point_num * sizeof(double));

//  Diagnostics
    cudaMalloc((void **)&Mh_d         , nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&W_d          , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Wh_d         , nvi* point_num *     sizeof(double));
    cudaMalloc((void **)&Rho_d        , nv * point_num   * sizeof(double));
    cudaMalloc((void **)&pressure_d   , nv * point_num   * sizeof(double));

//  Temperature
    cudaMalloc((void **)&temperature_d, nv * point_num *     sizeof(double));

//  Potential temperature
    cudaMalloc((void **)&pt_d         , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&pth_d        , nvi* point_num *     sizeof(double));

//  Entalphy
    cudaMalloc((void **)&h_d           , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&hh_d          , nvi * point_num *     sizeof(double));

//  Advection
    cudaMalloc((void **)&Adv_d        , nv * point_num * 3 * sizeof(double));

//  3D vector
    cudaMalloc((void **)&v_d          , nv * point_num * 3 * sizeof(double));

//  Effective gravity
    cudaMalloc((void **)&gtil_d        , nv * point_num * sizeof(double));
    cudaMalloc((void **)&gtilh_d       , nvi* point_num * sizeof(double));

//  Slow modes
    cudaMalloc((void **)&SlowMh_d        , nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&SlowWh_d        , nvi* point_num *     sizeof(double));
    cudaMalloc((void **)&SlowRho_d       , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Slowpressure_d  , nv * point_num *     sizeof(double));


//  Deviations
    cudaMalloc((void **)&pressures_d   , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Rhos_d        , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Mhs_d         , nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Ws_d          , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Whs_d         , nvi* point_num *     sizeof(double));



//  RK-Method
    cudaMalloc((void **)&pressurek_d   , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Rhok_d        , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Mhk_d         , nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&Wk_d          , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&Whk_d         , nvi* point_num *     sizeof(double));

//  Vertical integration
    cudaMalloc((void **)&Sp_d          , nv * point_num * sizeof(double));
    cudaMalloc((void **)&Sd_d          , nv * point_num * sizeof(double));

//  Diffusion
    cudaMalloc((void **)&Kdhz_d         ,nv *                 sizeof(double));
    cudaMalloc((void **)&Kdh4_d         ,nv *                 sizeof(double));
    cudaMalloc((void **)&DivM_d         ,nv * point_num * 3 * sizeof(double));
    cudaMalloc((void **)&diffpr_d       ,nv * point_num     * sizeof(double));
    cudaMalloc((void **)&diffmh_d       , 3 * nv * point_num* sizeof(double));
    cudaMalloc((void **)&diffw_d        , nv* point_num     * sizeof(double));
    cudaMalloc((void **)&diffrh_d       , nv * point_num    * sizeof(double));
    cudaMalloc((void **)&diff_d          , 6 * nv * point_num    * sizeof(double));
    cudaMalloc((void **)&divg_Mh_d       , 3 * nv * point_num    * sizeof(double));

//  Extras-nan
    cudaMalloc((void **)&check_d, sizeof (bool));

    cudaMalloc((void **)&vbar_d          , 3 * nv * point_num *sizeof(double));
    cudaMalloc((void **)&zonal_mean_tab_d, 2 * point_num * sizeof(int));
//  Rad Transfer
    cudaMalloc((void **)&fnet_up_d   , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&fnet_dn_d   , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&tau_d       , nv * point_num * 2 *  sizeof(double));

    cudaMalloc((void **)&phtemp      , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&thtemp      , nvi * point_num *     sizeof(double));
    cudaMalloc((void **)&ttemp       , nv * point_num *     sizeof(double));
    cudaMalloc((void **)&dtemp       , nv * point_num *     sizeof(double));
}

__host__ void ESP::InitialValues(bool rest          ,
                                 int glevel         ,
                                 double timestep_dyn,
                                 double A           ,
                                 double Cp          ,
                                 double P_Ref       ,
                                 double Gravit      ,
                                 double Omega       ,
                                 double Diffc       ,
                                 double kb          ,
                                 double Tmean       ,
                                 double Mmol        ,
                                 double mu          ,
                                 double Rd          ,
                                 int *zonal_mean_tab,
                                 bool sponge        ){
//
//  Description:
//
//  Set initial conditions.
//
//
//  Initial atmospheric conditions
    if(rest){
        for (int i = 0; i < point_num; i++ ){
//
//          Initial conditions for an isothermal Atmosphere
//
            double Ha = Rd * Tmean / Gravit;
            for (int lev = 0; lev < nv; lev++ ){
                pressure_h[i*nv + lev] = P_Ref*exp(-Altitude_h[lev] / Ha);
                temperature_h[i*nv + lev] = Tmean;
            }

            for (int lev = 0; lev < nv; lev++ ){
//              Density [kg/m3]
                Rho_h[i*nv + lev] = pressure_h[i*nv + lev] / (temperature_h[i*nv + lev] * Rd);

//              Momentum [kg/m3 m/s]
                Mh_h[i*3*nv + 3*lev + 0] = 0.0;
                Mh_h[i*3*nv + 3*lev + 1] = 0.0;
                Mh_h[i*3*nv + 3*lev + 2] = 0.0;

//              Vertical momentum [kg/m3 m/s]
                W_h[i*nv + lev] = 0.0;     // Center of the layer.
                Wh_h[i*(nv+1) + lev] = 0.0;// Layers interface.
            }
            Wh_h[i*(nv + 1) + nv] = 0.0;
        }
    }
    else{
//
//      Restart from an existing simulation.
        hid_t       file_id, dataset_id;

        char FILE_NAME1[160];

        sprintf(FILE_NAME1, "ifile/esp_initial.h5");
        file_id = H5Fopen(FILE_NAME1, H5F_ACC_RDONLY, H5P_DEFAULT);
//      Density
        dataset_id = H5Dopen(file_id, "/Rho",H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Rho_h);
//      Pressure
        dataset_id = H5Dopen(file_id, "/Pressure",H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure_h);
//      Horizontal momentum
        dataset_id = H5Dopen(file_id, "/Mh",H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Mh_h);
//      Vertical momentum
        dataset_id = H5Dopen(file_id, "/Wh",H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Wh_h);

        H5Dclose(dataset_id);
        H5Fclose(file_id);

        for(int lev = 0; lev < nv+1; lev++)
            for(int i = 0; i < point_num; i++)
                temperature_h[i*nv + lev] = pressure_h[i*nv + lev]/(Rd*Rho_h[i*nv + lev]);

        for(int i = 0; i < point_num; i++){
            for(int lev = 1; lev < nv; lev++){
                double xi  = Altitude_h[lev  ] ;
                double xim1= Altitudeh_h[lev ] ;
                double xip1= Altitudeh_h[lev +1  ] ;

                double a = (xi - xip1)/(xim1 -xip1);
                double b = (xi - xim1)/(xip1 -xim1);

                W_h[i*nv + lev] = Wh_h[i*(nv+1) + lev]*a + Wh_h[i*(nv+1) + lev+1]*b;
            }
        }
    }


//  Diffusion
//  Horizontal
    double *Kdhz_h, *Kdh4_h;
    Kdhz_h = new double[nv];
    Kdh4_h = new double[nv];
    for (int lev = 0; lev < nv; lev++ ){
//      Diffusion constant.
        double dbar = sqrt(2*M_PI/5)*A/(pow(2,glevel));
        Kdh4_h[lev] = Diffc*pow(dbar,4.)/timestep_dyn;
        Kdhz_h[lev] = Diffc*pow(dbar,4.)/timestep_dyn;
    }

//  Copy memory to the devide
    cudaMemcpy(point_local_d, point_local_h, 6 * point_num * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(maps_d, maps_h, (nl_region + 2)*(nl_region + 2)*nr * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(Altitude_d , Altitude_h , nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Altitudeh_d, Altitudeh_h, nvi* sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecoa_d  , nvecoa_h  , 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecti_d  , nvecti_h  , 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(nvecte_d  , nvecte_h  , 6 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(areasTr_d , areasTr_h , 6 *     point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(lonlat_d , lonlat_h , 2 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(func_r_d  , func_r_h  , 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(temperature_d, temperature_h, point_num * nv    *     sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Mh_d         , Mh_h         , point_num * nv    * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(W_d          , W_h          , point_num * nv    *     sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Wh_d         , Wh_h         , point_num * nvi   *     sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Rho_d        , Rho_h         , point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pressure_d   , pressure_h    , point_num * nv * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(div_d, div_h, 7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(grad_d,grad_h,7 * 3 * point_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdhz_d      ,Kdhz_h, nv     * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Kdh4_d      ,Kdh4_h, nv     * sizeof(double), cudaMemcpyHostToDevice);

//  Initialize arrays
    cudaMemset(Adv_d, 0, sizeof(double) * 3 * point_num * nv);
    cudaMemset(v_d  , 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(pt_d , 0, sizeof(double) * nv * point_num    );
    cudaMemset(pth_d, 0, sizeof(double) * nvi* point_num    );
    cudaMemset(SlowMh_d        , 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(SlowWh_d        , 0, sizeof(double) * nvi* point_num    );
    cudaMemset(SlowRho_d       , 0, sizeof(double) * nv * point_num    );
    cudaMemset(Slowpressure_d  , 0, sizeof(double) * nv * point_num    );
    cudaMemset(h_d        , 0, sizeof(double) * nv * point_num    );
    cudaMemset(hh_d       , 0, sizeof(double) * nvi * point_num    );
    cudaMemset(Rhos_d     , 0, sizeof(double) * nv * point_num    );
    cudaMemset(pressures_d, 0, sizeof(double) * nv * point_num    );
    cudaMemset(Mhs_d      , 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Ws_d       , 0, sizeof(double) * nv * point_num    );
    cudaMemset(Whs_d      , 0, sizeof(double) * nvi* point_num    );
    cudaMemset(gtil_d   , 0, sizeof(double) * nv * point_num);
    cudaMemset(gtilh_d  , 0, sizeof(double) * nvi* point_num);
    cudaMemset(Rhok_d     , 0, sizeof(double) * nv * point_num    );
    cudaMemset(pressurek_d, 0, sizeof(double) * nv * point_num    );
    cudaMemset(Mhk_d      , 0, sizeof(double) * nv * point_num * 3);
    cudaMemset(Wk_d       , 0, sizeof(double) * nv * point_num    );
    cudaMemset(Whk_d      , 0, sizeof(double) * nvi* point_num    );
    cudaMemset(Sp_d       , 0, sizeof(double) * point_num * nv);
    cudaMemset(Sd_d       , 0, sizeof(double) * point_num * nv);
    cudaMemset(DivM_d      , 0, sizeof(double) * point_num * 3 * nv);
    cudaMemset(diffpr_d    , 0, sizeof(double) * nv * point_num);
    cudaMemset(diffmh_d    , 0, sizeof(double) * 3 * nv * point_num);
    cudaMemset(diffw_d     , 0, sizeof(double) * nv * point_num);
    cudaMemset(diffrh_d    , 0, sizeof(double) * nv * point_num);
    cudaMemset(diff_d       , 0, sizeof(double) * 6 * nv * point_num);
    cudaMemset(divg_Mh_d    , 0, sizeof(double) * 3 * nv * point_num);

    if (sponge==true) cudaMemcpy(zonal_mean_tab_d      ,zonal_mean_tab, 2*point_num * sizeof(int), cudaMemcpyHostToDevice);

    delete [] Kdh4_h;
    delete [] Kdhz_h;

}

__host__ void ESP::RTSetup(double Tstar_           ,
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

__host__ ESP::~ESP(){

//
//  Description: Frees the memory space.
//
//  Host
    free(point_local_h);
    free(maps_h);
    free(lonlat_h);
    free(Altitude_h);
    free(Altitudeh_h);
    free(nvecoa_h);
    free(nvecti_h);
    free(nvecte_h);
    free(areasTr_h);
    free(div_h);
    free(grad_h);
    free(func_r_h);
    free(Rho_h);
    free(pressure_h);
    free(temperature_h);
    free(Mh_h);
    free(W_h);
    free(Wh_h);

//  Device
    cudaFree(point_local_d);
    cudaFree(maps_d);
    cudaFree(Altitude_d);
    cudaFree(Altitudeh_d);
    cudaFree(nvecoa_d);
    cudaFree(nvecti_d);
    cudaFree(nvecte_d);
    cudaFree(areasTr_d);
    cudaFree(lonlat_d);
    cudaFree(div_d);
    cudaFree(grad_d);
    cudaFree(func_r_d);
    cudaFree(Rho_d);
    cudaFree(pressure_d);
    cudaFree(temperature_d);
    cudaFree(W_d);
    cudaFree(Wh_d);
    cudaFree(h_d);
    cudaFree(hh_d);
    cudaFree(Adv_d);
    cudaFree(gtil_d);
    cudaFree(gtilh_d);
    cudaFree(v_d);
    cudaFree(pt_d);
    cudaFree(pth_d);
    cudaFree(SlowMh_d);
    cudaFree(SlowWh_d);
    cudaFree(SlowRho_d);
    cudaFree(Slowpressure_d);
    cudaFree(Rhok_d);
    cudaFree(pressurek_d);
    cudaFree(Mhk_d);
    cudaFree(Whk_d);
    cudaFree(Wk_d);
    cudaFree(Rhos_d);
    cudaFree(pressures_d);
    cudaFree(Mhs_d);
    cudaFree(Whs_d);
    cudaFree(Ws_d);
    cudaFree(Sd_d);
    cudaFree(Sp_d);
    cudaFree(Kdhz_d);
    cudaFree(Kdh4_d);
    cudaFree(DivM_d);
    cudaFree(diffpr_d);
    cudaFree(diffmh_d);
    cudaFree(diffw_d);
    cudaFree(diffrh_d);
    cudaFree(diff_d);
    cudaFree(divg_Mh_d);

    printf("\n\n Free memory!\n\n");
}
