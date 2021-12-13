// ==============================================================================
// This file is part of Alfrodull.
//
//     Alfrodull is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     Alfrodull is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     Alfrodull directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Method: Helios Two Stream algorithm
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// Code contributors: Urs Schroffenegger, Matej Malik
//
// History:
// Version Date       Comment
// ======= ====       =======
// 1.0     2020-07-15 First version
//
//
////////////////////////////////////////////////////////////////////////

#pragma once

#include <memory>

__device__ double analyt_planck(int n, double y1, double y2);


__device__ double power_int(double x, int i);

__host__ __device__ void thomas_solve(double4* A,
                                      double4* B,
                                      double4* C,
                                      double2* D,
                                      double4* C_prime,
                                      double2* D_prime,
                                      double2* X,
                                      int      N);

//***************************************************************************************************
__global__ void arrays_mean(double* array1, double* array2, double* array_out, int array_size);


//***************************************************************************************************

// first simple integration over weights
__global__ void integrate_val_band(double* val_wg,       // in
                                   double* val_band,     // out
                                   double* gauss_weight, // in
                                   int     nbin,
                                   int     num_val,
                                   int     ny);

// simple integration over bins/bands
__global__ void integrate_val_tot(double* val_tot,     // out
                                  double* val_band,    // in
                                  double* deltalambda, // in
                                  int     nbin,
                                  int     num_val);

std::shared_ptr<double[]>
integrate_band(double* val, double* gauss_weight, int num_val, int nbin, int ny);

std::shared_ptr<double[]> integrate_wg_band(double* val,
                                            double* gauss_weight,
                                            double* deltalambda,
                                            int     num_val,
                                            int     nbin,
                                            int     ny);
