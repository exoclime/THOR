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
// Description: Defines the icosahedral standard grid: variables and functions
//   
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

#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <math.h>

class Icogrid{

public:

    // Variables associated with the grid
    int    point_num;       
    int    nv       ;              
    int    nvi      ;             
    int    nl_region;
    int    nr       ;

    double *point_xyz ;  
    double *point_xyzq;

    int *pent_ind;
    int *halo    ;
    int *maps    ;

    double *Altitude; 
    double *Altitudeh; 

    double *lonlat;

    double *areas  ;
    double *areasT ;
    double *areasTr;

    double *nvec  ;
    double *nvecoa;
    double *nvecte;
    double *nvecti;

    double *div ;
    double *grad;

    int    *point_local;

    double *func_r;

    Icogrid (bool, double, int, int, double, double);
    
private:
    // Functions to build the grid
    void sphere_ico (double *, int, int, int, int, int, int, int *, int, int);                  
    void neighbors_indx ( int *, double *, int *, int);                                         
    void reorder_neighbors_indx (int *, double *, int*, int );                                  
    void neighbors_indx_pl(int *, double *, int *, int);                                        
    void reorder_neighbors_indx_pl(int *, double *, int*, int);                                 
    void generate_halos(int *, int *, int, int);
    void reorder_neighbors_indx_rhombi (int *, int *,int *, int, int, int, int, int);
    void produce_maps(int *, int *, int, int);
    void spring_dynamics(int *, int *, int, double, double *, int);
    void find_qpoints (int *, double *, double *, int *, int);
    void relocate_centres (int *, double *, double *, int *, int);
    void set_altitudes(double *, double *, double, int);
    void cart2sphe ( double *, double *, int); 
    void correct_xyz_points (double , double *, double *, int *, int);
    void control_areas(double *, double *, double *, int *, double *, double *, int *, int);
    void control_vec(double *, double *, double *, double *, double *, double *, double *, int *, int *, int);
    void compute_func(double *, double *, int);
    void div_operator(double *, double *, double *, double *, int *, int);
    void gra_operator(double *, double *, double *, double *, int *, int);
};

