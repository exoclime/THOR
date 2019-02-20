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
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////
#pragma once

#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER
#    define _USE_MATH_DEFINES
#endif

#include <math.h>


class Icogrid
{

public:
    // Variables associated with the grid
    int point_num; // number of horizontal points
    int nv;        // number of vertical layers
    int nvi;       // number of interfaces between leyers
    int nl_region; // number of  points in one side of a rhombus
    int nr;        // number of rombi

    int grid_level; // Number of recursive iterations to increase horizontal resolution.

    int nh; // number of points in halo

    double *point_xyz;
    double *point_xyzq;

    int *pent_ind;
    int *halo;
    int *maps; // Grid domains

    double *Altitude;  // Altitudes
    double *Altitudeh; // Altitude at the interfaces between layers

    double *lonlat; // Longitude and latitude of the grid points

    double *areas;
    double *areasT;  // Areas of the main cells
    double *areasTr; // Areas of the triangles

    double *nvec;
    double *nvecoa; // Normal vectors for diffusion 1
    double *nvecte; // Normal vectors for diffusion 3
    double *nvecti; // Normal vectors for diffusion 2

    double *mvec; // unit vectors used in curl operator

    double *div;   // Divergence operator
    double *grad;  // Gradient operator
    double *curlz; // vertical component of curl operator

    int *point_local; // First neighbours

    double *func_r; // Normalised vector

    int *zonal_mean_tab; //something something


    Icogrid(bool, double, int, int, int, double, double, bool, int *);
    void free_memory();


private:
    // Functions to build the grid
    void sphere_ico(double *, int, int, int, int, int, int, int *, int, int);
    void neighbors_indx(int *, double *, int *, int);
    void reorder_neighbors_indx(int *, double *, int *, int);
    void neighbors_indx_pl(int *, double *, int *, int);
    void reorder_neighbors_indx_pl(int *, double *, int *, int);
    void generate_halos(int *, int *, int, int);
    void reorder_neighbors_indx_rhombi(int *, int *, int *, int, int, int, int, int);
    void produce_maps(int *, int *, int, int);
    void spring_dynamics(int *, int *, int, double, double *, int);
    void find_qpoints(int *, double *, double *, int *, int);
    void relocate_centres(int *, double *, double *, int *, int);
    void set_altitudes(double *, double *, double, int);
    void cart2sphe(double *, double *, int);
    void correct_xyz_points(double, double *, double *, int *, int);
    void control_areas(double *, double *, double *, int *, double *, double *, int *, int);
    void control_vec(double *,
                     double *,
                     double *,
                     double *,
                     double *,
                     double *,
                     double *,
                     int *,
                     int *,
                     int,
                     double *);
    void compute_func(double *, double *, int);
    void div_operator(double *, double *, double *, double *, int *, int);
    void gra_operator(double *, double *, double *, double *, int *, int);
    void curlz_operator(double *, double *, double *, double *, int *, int);
    void zonal_mean_tab_f(int *, double *, int, int, int *);
};
