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
//
// Description: Sponge Layer
//
//
//
// Known limitations: None
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

__global__ void zonal_uv(double *M_d,
                         double *Rho_d,
                         int *   zonal_mean_tab_d,
                         double *lonlat_d,
                         int     num,
                         double *utmp,
                         double *vtmp,
                         int     max_count);

__global__ void
zonal_w(double *W_d, double *Rho_d, int *zonal_mean_tab_d, int num, double *wtmp, int max_count);

__global__ void zonal_temp(double *pressure_d,
                           double *Rho_d,
                           double *Tbar_d,
                           int *   zonal_mean_tab_d,
                           double *lonlat_d,
                           int     num,
                           double *Ttmp,
                           double *Rd_d,
                           int     max_count);

void print_vbar(double *vbar_h, int nlat, int nv);

__global__ void sponge_layer(double *M_d,
                             double *Rho_d,
                             double *W_d,
                             double *Wh_d,
                             double *pressure_d,
                             double *vbar_d,
                             double *Tbar_d,
                             int *   zonal_mean_tab_d,
                             double *lonlat_d,
                             double *Altitude_d,
                             double *Altitudeh_d,
                             double  Ruv,
                             double  Rw,
                             double  RT,
                             double  Rv_fac,
                             double  nsi,
                             bool    damp_uv_to_mean,
                             bool    damp_w_to_mean,
                             bool    implicit,
                             double  dt,
                             double *Rd_d,
                             int     nlat,
                             int     num,
                             int     nv,
                             bool    temp_sponge,
                             double *profx_dMh_d,
                             double *profx_dWh_d,
                             double *profx_dW_d,
                             double *profx_Qheat_d);
