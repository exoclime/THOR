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
//
//
// Description: Declares functions and main variables on the host and device
//   
// Method: - 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class ESP{

public:
///////////////////////////
//  General variable
    const int point_num   ;
    const int nv          ;
    const int nvi         ;
    const int nl_region   ;
    const int nr          ;
///////////////////////////
//  Host
    int *point_local_h    ;
    int *maps_h           ;

    double *lonlat_h      ;

    double *Altitude_h    ;
    double *Altitudeh_h   ;

    double *nvecoa_h      ;
    double *nvecti_h      ;
    double *nvecte_h      ;
    double *areasT_h      ;
    double *areasTr_h     ;
    
    double *div_h ;
    double *grad_h;

    double *func_r_h      ;

    double *Rho_h         ;
    double *pressure_h    ;
    double *temperature_h ;
    double *Mh_h          ;
    double *W_h           ;
    double *Wh_h          ;

    double *Kdhz_h        ;
    double *Kdh4_h        ;
    bool    check_h       ;
///////////////////////////
//  Device
    int *point_local_d    ;
    int *maps_d           ;

    double *Altitude_d    ;
    double *Altitudeh_d   ;

    double *nvecoa_d        ;
    double *nvecti_d      ;
    double *nvecte_d      ;
    double *areasTr_d     ;

    double *lonlat_d      ;

    double *div_d ;
    double *grad_d;

    double *func_r_d      ;

    double *temperature_d ;
    double *Mh_d          ;
    double *W_d           ;
    double *Wh_d          ;

    double *h_d           ;
    double *hh_d          ;

    double *Rho_d         ;
    double *pressure_d    ;

    double *Adv_d         ;

    double *SlowMh_d      ;
    double *SlowWh_d      ;
    double *SlowRho_d     ;
    double *Slowpressure_d;
    
    double *Rhos_d        ;
    double *pressures_d   ;
    double *Mhs_d         ;
    double *Ws_d          ;
    double *Whs_d         ;

    double *Rhok_d        ;
    double *pressurek_d   ;
    double *Mhk_d         ;
    double *Wk_d          ;
    double *Whk_d         ;

    double *v_d           ;
    double *pt_d          ;
    double *pth_d         ;
    
    double *gtil_d        ;
    double *gtilh_d       ;

    double *Sd_d          ;
    double *Sp_d          ;

    double *Kdhz_d        ;
    double *Kdh4_d        ;

    double *DivM_d        ;
    double *diffpr_d      ;
    double *diffmh_d      ;
    double *diffw_d       ;
    double *diffrh_d      ;
    
    double *diff_d        ;
    double *divg_Mh_d     ;
    bool   *check_d       ;
///////////////////////////

//  Functions
    ESP(int *   ,
        int *   ,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        double *,
        int     ,
        int     ,
        int     ,
        int     ,
        int     ); 

    void AllocData() ;

    void InitialValues(bool  , 
                       int   ,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double,
                       double);

    void Thor(double,
              bool  ,
              bool  ,
              double,
              double,
              double,
              double,
              double,
              double,
              double,
              double,
              double,
              bool  ,
              bool  );

    void ProfX(int   ,
               int   ,
               bool  ,
               double,
               double,
               double,
               double,
               double,
               double,
               double,
               double,
               double,
               double);

    void CopyToHost();

    void Output(int   ,
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                char* ,
                double);

    ~ESP();
};

