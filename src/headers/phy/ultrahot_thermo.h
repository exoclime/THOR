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
// Description: dry convective adjustment scheme
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

#include "directories.h"
#include <math.h>

#define Gibbs_Filename "src/headers/phy/GibbsH.txt"
#define RUNIV 8.3144621
#define massH 1.67372e-27
#define kBoltz 1.38064852e-23
#define Navo 6.0221409e23

// bool ESP::read_in_gibbs_H(int GibbsN);

__device__ __host__ double
           linear_interpolation(double x1, double x2, double y1, double y2, double x);

__device__ __host__ int locate_min_i(double *array_d, int N, double val);

__device__ __host__ int locate_max_i(double *array_d, int N, double val);

// double kprime(double *T, double *dG, double temperature, double pressure);

__device__ __host__ double
           chi_H_equilibrium(double *GibbsT, double *GibbsdG, int GibbsN, double temperature, double pressure);

__device__ __host__ double Rd_from_chi_H(double chi_H);

__device__ __host__ double heat_capacity_H2(double temperature);

__device__ __host__ double Cp_from_chi_H(double chi_H, double temperature);

__host__ double guillot_T(double pressure,
                          double mu,
                          double Teq,
                          double P_Ref,
                          double Gravit,
                          double Tint,
                          double f_lw,
                          double kappa_sw,
                          double kappa_lw);

__device__ double linear_interpolation_device(double x1, double x2, double y1, double y2, double x);

__device__ int locate_min_i_device(double *array_d, int N, double val);

__device__ int locate_max_i_device(double *array_d, int N, double val);

__device__ double chi_H_equilibrium_device(double *GibbsT,
                                           double *GibbsdG,
                                           int     GibbsN,
                                           double  temperature,
                                           double  pressure);

__device__ double Rd_from_chi_H_device(double chi_H);

__device__ double heat_capacity_H2_device(double temperature);

__device__ double Cp_from_chi_H_device(double chi_H, double temperature);


__global__ void update_temperature_Rd_Cp(double *temperature_d,
                                         double *Rd_d,
                                         double *Cp_d,
                                         double *pressure_d,
                                         double *Rho_d,
                                         double *GibbsT_d,
                                         double *GibbsdG_d,
                                         int     GibbsN,
                                         int     num,
                                         bool *  error);
