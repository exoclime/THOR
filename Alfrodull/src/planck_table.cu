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
// Planck table computation and interpolation
//
//
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


#include "planck_table.h"

#include "physics_constants.h"

#include "math_helpers.h"

// constructing a table with Planck function values for given wavelengths and in a suitable temperature range
__global__ void plancktable(double* planck_grid,
                            double* lambda_edge,
                            double* deltalambda,
                            int     nwave,
                            double  Tstar,
                            int     p_iter,
                            int     dim,
                            int     step) {
    // Wavenumber
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    // temperature
    int t = threadIdx.y + blockIdx.y * blockDim.y;

    if (x < nwave && t < (dim / 10 + 1)) {

        double T;
        double shifty;
        double D;
        double y_bot;
        double y_top;

        // building flexible temperature grid from '1 K' to 'dim * 2 - 1 K' at 'step K' resolution
        // and Tstar
        if (t < (dim / 10)) {
            T = (t + p_iter * (dim / 10)) * step + 1;
        }
        if (p_iter == 9) {
            if (t == dim / 10) {
                T = Tstar;
            }
        }

        planck_grid[x + (t + p_iter * (dim / 10)) * nwave] = 0.0;
        // analytical calculation, only for T > 0
        if (T > 0.01) {
            D = 2.0 * (power_int(KBOLTZMANN / HCONST, 3) * KBOLTZMANN * power_int(T, 4))
                / (CSPEED * CSPEED);
            y_top = HCONST * CSPEED / (lambda_edge[x + 1] * KBOLTZMANN * T);
            y_bot = HCONST * CSPEED / (lambda_edge[x] * KBOLTZMANN * T);

            // rearranging so that y_top < y_bot (i.e. wavelengths are always increasing)
            if (y_bot < y_top) {
                shifty = y_top;
                y_top  = y_bot;
                y_bot  = shifty;
            }

            for (int n = 1; n < 200; n++) {
                planck_grid[x + (t + p_iter * (dim / 10)) * nwave] +=
                    D * analyt_planck(n, y_bot, y_top);
            }
        }
        planck_grid[x + (t + p_iter * (dim / 10)) * nwave] /= deltalambda[x];
    }
}

// adjust the incident flux to correspond to the correct brightness temperature
__global__ void corr_inc_energy(double* planck_grid,
                                double* starflux,
                                double* deltalambda,
                                bool    realstar,
                                int     nwave,
                                double  Tstar,
                                int     dim) {

    int x = threadIdx.x + blockIdx.x * blockDim.x;

    if (x < nwave) {

        double num_flux = 0;

        if (realstar) {
            for (int xl = 0; xl < nwave; xl++) {

                num_flux += deltalambda[xl] * starflux[xl];
            }
        }
        else {
            for (int xl = 0; xl < nwave; xl++) {

                num_flux += deltalambda[xl] * PI * planck_grid[xl + dim * nwave];
            }
        }

        double theo_flux = STEFANBOLTZMANN * pow(Tstar, 4.0);

        double corr_factor = theo_flux / num_flux;
        if (x == 0) {
            if (corr_factor > 1)
                printf("\nEnergy budget corrected (increased) by %.2f percent.\n",
                       100.0 * (corr_factor - 1.0));
            if (corr_factor < 1)
                printf("\nEnergy budget corrected (decreased) by %.2f percent.\n",
                       100.0 * (1.0 - corr_factor));
        }
        if (realstar == 1) {

            starflux[x] *= corr_factor;
        }
        else {

            planck_grid[x + dim * nwave] *= corr_factor;
        }
    }
}


planck_table::planck_table() {
    // number of pre-tabulated temperature values for the planck table
    dim = 8000;
    // temperature step for the planck table. e.g. dim = 10000 and step = 2 will give a table from 1 K to 19999 K in 2 K steps
    step = 2;
}

planck_table::~planck_table() {
}

void planck_table::construct_planck_table(
    double* lambda_edge, // linked to opacity tables binning: opacities.dev_opac_interwave
    double* deltalambda, // linked to opacity tables binning: opacities.dev_opac_deltawave
    int     nwave,       // linked to opacity tables binning: nbin
    double  Tstar_) {
    nplanck_grid = (dim + 1) * nwave;
    Tstar        = Tstar_;

    planck_grid.allocate(nplanck_grid);

    dim3 grid((int(nwave) + 15) / 16, (int(dim / 10.0 + 1.0) + 15) / 16, 1);
    dim3 block(16, 16, 1);

    for (int p_iter = 0; p_iter < 10; p_iter++) {
        plancktable<<<grid, block>>>(
            *planck_grid, lambda_edge, deltalambda, nwave, Tstar, p_iter, dim, step);
        cudaDeviceSynchronize();
    }
}
