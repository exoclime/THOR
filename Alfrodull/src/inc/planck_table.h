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

#include "cuda_device_memory.h"

// adjust the incident flux to correspond to the correct brightness temperature
__global__ void corr_inc_energy(double* planck_grid,
                                double* starflux,
                                double* deltalambda,
                                bool    realstar,
                                int     nwave,
                                double  Tstar,
                                int     dim);

class planck_table
{
public:
    planck_table();

    ~planck_table();


    void construct_planck_table(double* lambda_edge, double* deltalambda, int nwave, double Tstar);

    cuda_device_memory<double> planck_grid;

    int    dim   = 0;
    int    step  = 0;
    double Tstar = 0;
    int    nplanck_grid;
};
