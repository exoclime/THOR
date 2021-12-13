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
// kernels doing the flux integration, beam computation, iterative and thomas solver
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


#include "calculate_physics.h"
#include "physics_constants.h"

#include "math_helpers.h"

#include <stdio.h>

// calculates the integrated upwards and downwards fluxes
/*
  NOTE: called as this:

    dim3 threadsPerBlock(1, 1, 1); <- grid dim, defines BlockIdx
    dim3 numBlocks(32, 4, 8);      <- threads in a block, -> defines blockDim and threadIdx

    int nbin = opacities.nbin;
    int ny   = opacities.ny;

    //printf("Running Alfrodull Wrapper for integrate flux\n");
    integrate_flux_double<<<threadsPerBlock, numBlocks>>>(deltalambda,
*/
/*
// Old version, not used, not adapted to multi-column version.
// is slow due to usage of blocking AtomicAdd which effectively serialises
// additions.
__global__ void integrate_flux_double(double* deltalambda,  // in
                                      double* F_down_tot,   // out
                                      double* F_up_tot,     // out
                                      double* F_dir_tot,    // out
                                      double* F_net,        // out
                                      double* F_down_wg,    // in
                                      double* F_up_wg,      // in
                                      double* F_dir_wg,     // in
                                      double* F_down_band,  // out
                                      double* F_up_band,    // out
                                      double* F_dir_band,   // out
                                      double* gauss_weight, // in
                                      int     nbin,
                                      int     numinterfaces,
                                      int     ny) {

    int x = threadIdx.x;
    int y = threadIdx.y;
    int i = threadIdx.z;
    // set memory to 0.

    if (y == 0) {
        while (i < numinterfaces) {
            while (x < nbin) {

                F_up_tot[i]   = 0;
                F_down_tot[i] = 0;

                F_dir_band[x + nbin * i]  = 0;
                F_up_band[x + nbin * i]   = 0;
                F_down_band[x + nbin * i] = 0;

                x += blockDim.x;
            }
            x = threadIdx.x;
            i += blockDim.z;
        }
    }
    __syncthreads();

    i = threadIdx.z;

    while (i < numinterfaces) {
        while (y < ny) {
            while (x < nbin) {

                atomicAdd_double(&(F_dir_band[x + nbin * i]),
                                 0.5 * gauss_weight[y] * F_dir_wg[y + ny * x + ny * nbin * i]);
                atomicAdd_double(&(F_up_band[x + nbin * i]),
                                 0.5 * gauss_weight[y] * F_up_wg[y + ny * x + ny * nbin * i]);
                atomicAdd_double(&(F_down_band[x + nbin * i]),
                                 0.5 * gauss_weight[y] * F_down_wg[y + ny * x + ny * nbin * i]);

                x += blockDim.x;
            }
            x = threadIdx.x;
            y += blockDim.y;
        }
        y = threadIdx.y;
        i += blockDim.z;
    }
    __syncthreads();

    i = threadIdx.z;

    if (y == 0) {
        while (i < numinterfaces) {
            while (x < nbin) {

                atomicAdd_double(&(F_up_tot[i]), F_up_band[x + nbin * i] * deltalambda[x]);
                atomicAdd_double(&(F_down_tot[i]),
                                 (F_dir_band[x + nbin * i] + F_down_band[x + nbin * i])
                                     * deltalambda[x]);

                atomicAdd_double(&(F_dir_tot[i]), F_dir_band[x + nbin * i] * deltalambda[x]);

                x += blockDim.x;
            }
            x = threadIdx.x;
            i += blockDim.z;
        }
    }
    __syncthreads();

    i = threadIdx.z;

    if (x == 0 && y == 0) {
        while (i < numinterfaces) {
            F_net[i] = F_up_tot[i] - F_down_tot[i];
            i += blockDim.z;
        }
    }
}
*/
/*
numinterfaces: levels
nbin:          frequency bins
ny:            weights in frequency bins

must sum:
* over weighted ny into frequency bins (down_band, up_band. dir_band)
* over frequency bins into total flux per levels (up_tot, down_tot, net)
*/
// first simple integration over weights
__global__ void integrate_flux_band(double* F_down_wg,         // in
                                    double* F_up_wg,           // in
                                    double* F_dir_wg,          // in
                                    double* F_down_band,       // out
                                    double* F_up_band,         // out
                                    double* F_dir_band,        // out
                                    double* F_up_TOA_spectrum, // out
                                    double* gauss_weight,      // in
                                    int     nbin,
                                    int     numinterfaces,
                                    int     ny) {

    int interface_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int bin_idx       = blockIdx.y * blockDim.y + threadIdx.y;
    int c             = blockIdx.z;


    if (interface_idx < numinterfaces && bin_idx < nbin) {
        // set memory to 0.
        F_dir_band[bin_idx + nbin * interface_idx + c * nbin * numinterfaces] = 0;

        F_up_band[bin_idx + nbin * interface_idx + c * nbin * numinterfaces]   = 0;
        F_down_band[bin_idx + nbin * interface_idx + c * nbin * numinterfaces] = 0;

        int bin_offset = bin_idx + nbin * interface_idx + c * nbin * numinterfaces;

        for (int y = 0; y < ny; y++) {
            double w = gauss_weight[y];
            int    weight_offset =
                y + ny * bin_idx + ny * nbin * interface_idx + c * nbin * ny * numinterfaces;

            F_dir_band[bin_offset] += 0.5 * w * F_dir_wg[weight_offset];
            F_up_band[bin_offset] += 0.5 * w * F_up_wg[weight_offset];
            F_down_band[bin_offset] += 0.5 * w * F_down_wg[weight_offset];
        }

        // copy TOA spectrum to spectra array for TOA at each col.
        if (interface_idx == numinterfaces - 1)
            F_up_TOA_spectrum[bin_idx + c * nbin] = F_up_band[bin_offset];
    }
}

// simple integration over bins/bands
__global__ void integrate_flux_tot(double* deltalambda, // in
                                   double* F_down_tot,  // out
                                   double* F_up_tot,    // out
                                   double* F_dir_tot,   // out
                                   double* F_net,       // out
                                   double* F_down_band, // out
                                   double* F_up_band,   // out
                                   double* F_dir_band,  // out
                                   int     nbin,
                                   int     numinterfaces) {


    int interface_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int c             = blockIdx.z;

    if (interface_idx < numinterfaces) {

        F_up_tot[interface_idx + c * numinterfaces]   = 0;
        F_down_tot[interface_idx + c * numinterfaces] = 0;
        F_dir_tot[interface_idx + c * numinterfaces]  = 0;

        for (int bin = 0; bin < nbin; bin++) {
            int band_idx = interface_idx * nbin + bin + c * numinterfaces * nbin;
            F_up_tot[interface_idx + c * numinterfaces] += F_up_band[band_idx] * deltalambda[bin];
            F_down_tot[interface_idx + c * numinterfaces] +=
                (F_down_band[band_idx] + F_dir_band[band_idx]) * deltalambda[bin];
            F_dir_tot[interface_idx + c * numinterfaces] += F_dir_band[band_idx] * deltalambda[bin];
        }

        __syncthreads();
        F_net[interface_idx + c * numinterfaces] = F_up_tot[interface_idx + c * numinterfaces]
                                                   - F_down_tot[interface_idx + c * numinterfaces];
    }
}


// calculates the direct beam flux with geometric zenith angle correction, isothermal version
__global__ void fdir_iso(double* F_dir_wg,       // out
                         double* planckband_lay, // in
                         double* delta_tau_wg,   // in
                         double* z_lay,          // in
                         double* mu_star_cols,
                         double  mu_star_limit,
                         double  R_planet,
                         double  R_star,
                         double  a,
                         bool    dir_beam,
                         bool    geom_zenith_corr,
                         int     ninterface,
                         int     nbin,
                         int     ny,
                         int     num_cols) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int x = threadIdx.y + blockIdx.y * blockDim.y;

    // compound column and weight index
    int z = threadIdx.z + blockIdx.z * blockDim.z;
    int y = z % ny;
    // column
    int c      = z / ny;
    int nlayer = ninterface - 1;

    if (i < ninterface && x < nbin && y < ny && c < num_cols) {
        double mu_star = mu_star_cols[c];

        // the stellar intensity at TOA
        double I_dir =
            ((R_star / a) * (R_star / a)) * PI
            * planckband_lay[(ninterface - 1) + x * (ninterface - 1 + 2) + c * nbin * (nlayer + 2)];

        double f_out = 0.0;
        double mu_star_layer_j;

        // initialize each flux value
        if (dir_beam) {
            if (mu_star > mu_star_limit)
                F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = 0.0;
            else {
                f_out = -mu_star * I_dir;

                for (int j = ninterface - 2; j >= i; j--) {
                    // flux values lower that TOA will now be attenuated depending on their location

                    if (geom_zenith_corr) {
                        mu_star_layer_j =
                            -sqrt(1.0
                                  - pow((R_planet + z_lay[i]) / (R_planet + z_lay[j]), 2.0)
                                        * (1.0 - pow(mu_star, 2.0)));
                    }
                    else {
                        mu_star_layer_j = mu_star;
                    }

                    // direct stellar flux
                    f_out *= exp(delta_tau_wg[y + ny * x + ny * nbin * j + c * nlayer * ny * nbin]
                                 / mu_star_layer_j);
                }
                F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = f_out;
            }
        }
        else
            F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = 0.0;
    }
}


// calculates the direct beam flux with geometric zenith angle correction, non-isothermal version
__global__ void fdir_noniso(double* F_dir_wg,
                            double* Fc_dir_wg,
                            double* planckband_lay,
                            double* delta_tau_wg_upper,
                            double* delta_tau_wg_lower,
                            double* z_lay,
                            double* mu_star_cols,
                            double  mu_star_limit,
                            double  R_planet,
                            double  R_star,
                            double  a,
                            bool    dir_beam,
                            bool    geom_zenith_corr,
                            int     ninterface,
                            int     nbin,
                            int     ny,
                            int     num_cols) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int x = threadIdx.y + blockIdx.y * blockDim.y;

    // compound column and weight index
    int z = threadIdx.z + blockIdx.z * blockDim.z;
    int y = z % ny;
    // column
    int c      = z / ny;
    int nlayer = ninterface - 1;

    if (i < ninterface && x < nbin && y < ny && c < num_cols) {
        double mu_star = mu_star_cols[c];

        // the stellar intensity at TOA
        double I_dir =
            ((R_star / a) * (R_star / a)) * PI
            * planckband_lay[(ninterface - 1) + x * (ninterface - 1 + 2) + c * nbin * (nlayer + 2)];
        double f_out  = 0.0;
        double fc_out = 0.0;
        double mu_star_layer_j;

        // initialize each flux value
        if (dir_beam) {
            if (mu_star > mu_star_limit) {
                Fc_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = 0.0;
                F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin]  = 0.0;
            }
            else {
                f_out = -mu_star * I_dir;

                // flux values lower that TOA will now be attenuated depending on their location
                for (int j = ninterface - 2; j >= i; j--) {

                    if (geom_zenith_corr) {
                        mu_star_layer_j =
                            -sqrt(1.0
                                  - pow((R_planet + z_lay[i]) / (R_planet + z_lay[j]), 2.0)
                                        * (1.0 - pow(mu_star, 2.0)));
                    }
                    else {
                        mu_star_layer_j = mu_star;
                    }

                    double delta_tau =
                        delta_tau_wg_upper[y + ny * x + ny * nbin * j + c * nlayer * ny * nbin]
                        + delta_tau_wg_lower[y + ny * x + ny * nbin * j + c * nlayer * ny * nbin];

                    // direct stellar flux
                    fc_out = f_out
                             * exp(delta_tau_wg_upper[y + ny * x + ny * nbin * j
                                                      + c * nlayer * ny * nbin]
                                   / mu_star_layer_j);
                    f_out *= exp(delta_tau / mu_star_layer_j);
                }

                Fc_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = fc_out;
                F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin]  = f_out;
            }
        }
        else {
            Fc_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin] = 0.0;
            F_dir_wg[y + ny * x + ny * nbin * i + c * ninterface * ny * nbin]  = 0.0;
        }
    }
}


// calculation of the spectral fluxes, isothermal case with emphasis on on-the-fly calculations
// HELIOS iterative solver version
__global__ void fband_iso_notabu(double* F_down_wg_,      // out
                                 double* F_up_wg_,        // out
                                 double* F_dir_wg_,       // in
                                 double* planckband_lay_, // in
                                 double* w_0_,            // in
                                 double* M_term_,         // in
                                 double* N_term_,         // in
                                 double* P_term_,         // in
                                 double* G_plus_,         // in
                                 double* G_minus_,        // in
                                 double* g_0_tot_,        // in (clouds)
                                 double* surface_albedo,
                                 bool    singlewalk,
                                 double  Rstar,
                                 double  a,
                                 int     numinterfaces,
                                 int     nbin,
                                 double  f_factor,
                                 double* mu_star_cols,
                                 int     ny,
                                 int     num_cols,
                                 double  epsi,
                                 bool    dir_beam,
                                 bool    clouds,
                                 bool    scat_corr,
                                 bool    debug,
                                 double  i2s_transition) {

    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;

    int c = blockIdx.z;
    if (x < nbin && y < ny && c < num_cols) {
        int nlayer = numinterfaces - 1;
        // Start by applying column offset
        double* F_down_wg      = &(F_down_wg_[c * numinterfaces * nbin * ny]);
        double* F_up_wg        = &(F_up_wg_[c * numinterfaces * nbin * ny]);
        double* F_dir_wg       = &(F_dir_wg_[c * numinterfaces * nbin * ny]);
        double* planckband_lay = &(planckband_lay_[c * (nlayer + 2) * nbin]);
        double* w_0            = &(w_0_[c * nlayer * nbin * ny]);
        double* M_term         = &(M_term_[c * nlayer * nbin * ny]);
        double* N_term         = &(N_term_[c * nlayer * nbin * ny]);
        double* P_term         = &(P_term_[c * nlayer * nbin * ny]);
        double* G_plus         = &(G_plus_[c * nlayer * nbin * ny]);
        double* G_minus        = &(G_minus_[c * nlayer * nbin * ny]);
        double* g_0_tot        = &(g_0_tot_[c * nlayer * nbin * ny]);


        double mu_star = mu_star_cols[c];
        double w0;
        double M;
        double N;
        double P;
        double G_pl;
        double G_min;
        double g0;

        double E;

        double flux_terms;
        double planck_terms;
        double direct_terms;

        // calculation of downward fluxes from TOA to BOA
        for (int i = numinterfaces - 1; i >= 0; i--) {

            // TOA boundary -- incoming stellar flux
            if (i == numinterfaces - 1) {
                if (dir_beam)
                    F_down_wg[y + ny * x + ny * nbin * i] = 0.0;
                else
                    F_down_wg[y + ny * x + ny * nbin * i] =
                        f_factor * ((Rstar / a) * (Rstar / a)) * PI
                        * planckband_lay[i + x * (numinterfaces - 1 + 2)];
            }
            else {
                w0    = w_0[y + ny * x + ny * nbin * i];
                M     = M_term[y + ny * x + ny * nbin * i];
                N     = N_term[y + ny * x + ny * nbin * i];
                P     = P_term[y + ny * x + ny * nbin * i];
                G_pl  = G_plus[y + ny * x + ny * nbin * i];
                G_min = G_minus[y + ny * x + ny * nbin * i];
                g0    = g_0_tot[y + ny * x + ny * nbin * i];

                // DBG:
                // printf("%g %g %g %g %g %g %g\n",
                //        w0,
                //        M,
                //        N,
                //        P,
                //        G_pl,
                //        G_min,
                //        g0);
                // improved scattering correction factor E
                E = 1.0;
                if (scat_corr) {
                    E = E_parameter(w0, g0, i2s_transition);
                }

                // isothermal solution
                flux_terms = P * F_down_wg[y + ny * x + ny * nbin * (i + 1)]
                             - N * F_up_wg[y + ny * x + ny * nbin * i];

                planck_terms = planckband_lay[i + x * (numinterfaces - 1 + 2)] * (N + M - P);

                direct_terms = -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min * M + G_pl * N)
                               + F_dir_wg[y + ny * x + ny * nbin * (i + 1)] * P * G_min;

                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                F_down_wg[y + ny * x + ny * nbin * i] =
                    1.0 / M
                    * (flux_terms + 2.0 * PI * epsi * (1.0 - w0) / (E - w0) * planck_terms
                       + direct_terms);

                // DBG:
                //		if (isnan(F_down_wg[y + ny * x + ny * nbin * i]))
                // printf("NaN: %d %d %d - Fdip1: %g Fui: %g ft: %g pl: %g dt: %g\n",
                //        y,
                //        x,
                //        i,
                //        F_down_wg[y + ny * x + ny * nbin * (i + 1)],
                //        F_up_wg[y + ny * x + ny * nbin * i],
                //        flux_terms,
                //        planck_terms,
                //        direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (F_down_wg[y + ny * x + ny * nbin * i] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i,
                               x,
                               y);
                }
            }
        }

        // calculation of upward fluxes from BOA to TOA
        for (int i = 0; i < numinterfaces; i++) {

            // BOA boundary -- surface emission and reflection
            if (i == 0) {
                double w0_N = w_0[y + ny * x + ny * nbin * 0];
                double g0_N = g_0_tot[y + ny * x + ny * nbin * 0];

                // improved scattering correction factor E
                double E_N = 1.0;
                if (scat_corr) {
                    E_N = E_parameter(w0_N, g0_N, i2s_transition);
                }

                double BOA_part =
                    (1 - surface_albedo[x]) * PI * (1.0 - w0_N) / (E_N - w0_N)
                    * planckband_lay[numinterfaces
                                     + x * (numinterfaces - 1 + 2)]; // ghost layer plank emission
                double reflected_part = surface_albedo[x]
                                        * (F_dir_wg[y + ny * x + ny * nbin * i]
                                           + F_down_wg[y + ny * x + ny * nbin * i]);

                F_up_wg[y + ny * x + ny * nbin * i] = BOA_part + reflected_part;
            }
            else {
                w0    = w_0[y + ny * x + ny * nbin * (i - 1)];
                M     = M_term[y + ny * x + ny * nbin * (i - 1)];
                N     = N_term[y + ny * x + ny * nbin * (i - 1)];
                P     = P_term[y + ny * x + ny * nbin * (i - 1)];
                G_pl  = G_plus[y + ny * x + ny * nbin * (i - 1)];
                G_min = G_minus[y + ny * x + ny * nbin * (i - 1)];
                g0    = g_0_tot[y + ny * x + ny * nbin * (i - 1)];

                // improved scattering correction factor E
                E = 1.0;
                if (scat_corr) {
                    E = E_parameter(w0, g0, i2s_transition);
                }

                // isothermal solution
                flux_terms = P * F_up_wg[y + ny * x + ny * nbin * (i - 1)]
                             - N * F_down_wg[y + ny * x + ny * nbin * i];

                planck_terms = planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * (N + M - P);

                direct_terms = -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min * N + G_pl * M)
                               + F_dir_wg[y + ny * x + ny * nbin * (i - 1)] * P * G_pl;
                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                F_up_wg[y + ny * x + ny * nbin * i] =
                    1.0 / M
                    * (flux_terms + 2.0 * PI * epsi * (1.0 - w0) / (E - w0) * planck_terms
                       + direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (F_up_wg[y + ny * x + ny * nbin * i] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i,
                               x,
                               y);
                }
            }
        }
    }
}

// calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
// HELIOS iterative solver version
__global__ void fband_noniso_notabu(double* F_down_wg_,
                                    double* F_up_wg_,
                                    double* Fc_down_wg_,
                                    double* Fc_up_wg_,
                                    double* F_dir_wg_,
                                    double* Fc_dir_wg_,
                                    double* planckband_lay_,
                                    double* planckband_int_,
                                    double* w_0_upper_,
                                    double* w_0_lower_,
                                    double* delta_tau_wg_upper_,
                                    double* delta_tau_wg_lower_,
                                    double* M_upper_,
                                    double* M_lower_,
                                    double* N_upper_,
                                    double* N_lower_,
                                    double* P_upper_,
                                    double* P_lower_,
                                    double* G_plus_upper_,
                                    double* G_plus_lower_,
                                    double* G_minus_upper_,
                                    double* G_minus_lower_,
                                    double* g_0_tot_upper_,
                                    double* g_0_tot_lower_,
                                    double* surface_albedo,
                                    bool    singlewalk,
                                    double  Rstar,
                                    double  a,
                                    int     numinterfaces,
                                    int     nbin,
                                    double  f_factor,
                                    double* mu_star_cols,
                                    int     ny,
                                    int     num_cols,
                                    double  epsi,
                                    double  delta_tau_limit,
                                    bool    dir_beam,
                                    bool    clouds,
                                    bool    scat_corr,
                                    bool    debug,
                                    double  i2s_transition) {

    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int c = blockIdx.z;

    int nlayer = numinterfaces - 1;

    if (x < nbin && y < ny && c < num_cols) {
        double* F_down_wg          = &(F_down_wg_[c * numinterfaces * nbin * ny]);
        double* F_up_wg            = &(F_up_wg_[c * numinterfaces * nbin * ny]);
        double* Fc_down_wg         = &(Fc_down_wg_[c * numinterfaces * nbin * ny]);
        double* Fc_up_wg           = &(Fc_up_wg_[c * numinterfaces * nbin * ny]);
        double* F_dir_wg           = &(F_dir_wg_[c * numinterfaces * nbin * ny]);
        double* Fc_dir_wg          = &(Fc_dir_wg_[c * numinterfaces * nbin * ny]);
        double* planckband_lay     = &(planckband_lay_[c * (nlayer + 2) * nbin]);
        double* planckband_int     = &(planckband_int_[c * (nlayer + 1) * nbin]);
        double* w_0_upper          = &(w_0_upper_[c * nlayer * nbin * ny]);
        double* w_0_lower          = &(w_0_lower_[c * nlayer * nbin * ny]);
        double* delta_tau_wg_upper = &(delta_tau_wg_upper_[c * nlayer * nbin * ny]);
        double* delta_tau_wg_lower = &(delta_tau_wg_lower_[c * nlayer * nbin * ny]);
        double* M_upper            = &(M_upper_[c * nlayer * nbin * ny]);
        double* M_lower            = &(M_lower_[c * nlayer * nbin * ny]);
        double* N_upper            = &(N_upper_[c * nlayer * nbin * ny]);
        double* N_lower            = &(N_lower_[c * nlayer * nbin * ny]);
        double* P_upper            = &(P_upper_[c * nlayer * nbin * ny]);
        double* P_lower            = &(P_lower_[c * nlayer * nbin * ny]);
        double* G_plus_upper       = &(G_plus_upper_[c * nlayer * nbin * ny]);
        double* G_plus_lower       = &(G_plus_lower_[c * nlayer * nbin * ny]);
        double* G_minus_upper      = &(G_minus_upper_[c * nlayer * nbin * ny]);
        double* G_minus_lower      = &(G_minus_lower_[c * nlayer * nbin * ny]);
        double* g_0_tot_upper      = &(g_0_tot_upper_[c * nlayer * nbin * ny]);
        double* g_0_tot_lower      = &(g_0_tot_lower_[c * nlayer * nbin * ny]);

        double mu_star = mu_star_cols[c];

        double w0_up;
        double del_tau_up;
        double M_up;
        double N_up;
        double P_up;
        double G_pl_up;
        double G_min_up;
        double g0_up;
        double E_up;

        double w0_low;
        double del_tau_low;
        double M_low;
        double N_low;
        double P_low;
        double G_pl_low;
        double G_min_low;
        double g0_low;
        double E_low;

        double flux_terms;
        double planck_terms;
        double direct_terms;

        // calculation of downward fluxes from TOA to BOA
        for (int i = numinterfaces - 1; i >= 0; i--) {

            // TOA boundary -- incoming stellar flux
            if (i == numinterfaces - 1) {
                if (dir_beam)
                    F_down_wg[y + ny * x + ny * nbin * i] = 0.0;
                else
                    F_down_wg[y + ny * x + ny * nbin * i] =
                        f_factor * ((Rstar / a) * (Rstar / a)) * PI
                        * planckband_lay[i + x * (numinterfaces - 1 + 2)];
            }
            else {
                // upper part of layer quantities
                w0_up      = w_0_upper[y + ny * x + ny * nbin * i];
                del_tau_up = delta_tau_wg_upper[y + ny * x + ny * nbin * i];
                M_up       = M_upper[y + ny * x + ny * nbin * i];
                N_up       = N_upper[y + ny * x + ny * nbin * i];
                P_up       = P_upper[y + ny * x + ny * nbin * i];
                G_pl_up    = G_plus_upper[y + ny * x + ny * nbin * i];
                G_min_up   = G_minus_upper[y + ny * x + ny * nbin * i];
                g0_up      = g_0_tot_upper[y + ny * x + ny * nbin * i];

                // lower part of layer quantities
                w0_low      = w_0_lower[y + ny * x + ny * nbin * i];
                del_tau_low = delta_tau_wg_lower[y + ny * x + ny * nbin * i];
                M_low       = M_lower[y + ny * x + ny * nbin * i];
                N_low       = N_lower[y + ny * x + ny * nbin * i];
                P_low       = P_lower[y + ny * x + ny * nbin * i];
                G_pl_low    = G_plus_lower[y + ny * x + ny * nbin * i];
                G_min_low   = G_minus_lower[y + ny * x + ny * nbin * i];
                g0_low      = g_0_tot_lower[y + ny * x + ny * nbin * i];

                // improved scattering correction factor E
                E_up  = 1.0;
                E_low = 1.0;

                // improved scattering correction disabled for the following terms -- at least for the moment
                if (scat_corr) {
                    E_up  = E_parameter(w0_up, g0_up, i2s_transition);
                    E_low = E_parameter(w0_low, g0_low, i2s_transition);
                }

                // upper part of layer calculations
                if (del_tau_up < delta_tau_limit) {
                    // the isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = (planckband_int[(i + 1) + x * numinterfaces]
                                    + planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                   / 2.0 * (N_up + M_up - P_up);
                }
                else {
                    // the non-isothermal solution -- standard case
                    double pgrad_up = (planckband_lay[i + x * (numinterfaces - 1 + 2)]
                                       - planckband_int[(i + 1) + x * numinterfaces])
                                      / del_tau_up;

                    planck_terms =
                        planckband_lay[i + x * (numinterfaces - 1 + 2)] * (M_up + N_up)
                        - planckband_int[(i + 1) + x * numinterfaces] * P_up
                        + epsi / (E_up * (1.0 - w0_up * g0_up)) * (P_up - M_up + N_up) * pgrad_up;
                }
                flux_terms = P_up * F_down_wg[y + ny * x + ny * nbin * (i + 1)]
                             - N_up * Fc_up_wg[y + ny * x + ny * nbin * i];

                direct_terms =
                    -Fc_dir_wg[y + ny * x + ny * nbin * i] * (G_min_up * M_up + G_pl_up * N_up)
                    + F_dir_wg[y + ny * x + ny * nbin * (i + 1)] * G_min_up * P_up;
                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                Fc_down_wg[y + ny * x + ny * nbin * i] =
                    1.0 / M_up
                    * (flux_terms + 2.0 * PI * epsi * (1.0 - w0_up) / (E_up - w0_up) * planck_terms
                       + direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (Fc_down_wg[y + ny * x + ny * nbin * i] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i,
                               x,
                               y);
                }

                // lower part of layer calculations
                if (del_tau_low < delta_tau_limit) {
                    // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = (planckband_int[i + x * numinterfaces]
                                    + planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                   / 2.0 * (N_low + M_low - P_low);
                }
                else {
                    // non-isothermal solution -- standard case
                    double pgrad_low = (planckband_int[i + x * numinterfaces]
                                        - planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                       / del_tau_low;

                    planck_terms = planckband_int[i + x * numinterfaces] * (M_low + N_low)
                                   - planckband_lay[i + x * (numinterfaces - 1 + 2)] * P_low
                                   + epsi / (E_low * (1.0 - w0_low * g0_low))
                                         * (P_low - M_low + N_low) * pgrad_low;
                }
                flux_terms = P_low * Fc_down_wg[y + ny * x + ny * nbin * i]
                             - N_low * F_up_wg[y + ny * x + ny * nbin * i];

                direct_terms =
                    -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min_low * M_low + G_pl_low * N_low)
                    + Fc_dir_wg[y + ny * x + ny * nbin * i] * P_low * G_min_low;
                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                F_down_wg[y + ny * x + ny * nbin * i] =
                    1.0 / M_low
                    * (flux_terms
                       + 2.0 * PI * epsi * (1.0 - w0_low) / (E_low - w0_low) * planck_terms
                       + direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (F_down_wg[y + ny * x + ny * nbin * i] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i,
                               x,
                               y);
                }
            }
        }

        __syncthreads();

        // calculation of upward fluxes from BOA to TOA
        for (int i = 0; i < numinterfaces; i++) {

            // BOA boundary -- surface emission and reflection
            if (i == 0) {
                double reflected_part = surface_albedo[x]
                                        * (F_dir_wg[y + ny * x + ny * nbin * i]
                                           + F_down_wg[y + ny * x + ny * nbin * i]);

                double w0_N = w_0_lower[y + ny * x + ny * nbin * 0];
                double g0_N = g_0_tot_lower[y + ny * x + ny * nbin * 0];

                // improved scattering correction factor E
                double E_N = 1.0;
                if (scat_corr) {
                    E_N = E_parameter(w0_N, g0_N, i2s_transition);
                }

                // this is the surface/BOA emission. it correctly includes the emissivity e = (1 - albedo)
                double BOA_part = (1 - surface_albedo[x]) * PI * (1.0 - w0_N) / (E_N - w0_N)
                                  * planckband_lay[numinterfaces + x * (numinterfaces - 1 + 2)];
                // double BOA_part = 0.0;

                F_up_wg[y + ny * x + ny * nbin * i] =
                    reflected_part
                    + BOA_part; // internal_part consists of the internal heat flux plus the surface/BOA emission
            }
            else {
                // lower part of layer quantities
                w0_low      = w_0_lower[y + ny * x + ny * nbin * (i - 1)];
                del_tau_low = delta_tau_wg_lower[y + ny * x + ny * nbin * (i - 1)];
                M_low       = M_lower[y + ny * x + ny * nbin * (i - 1)];
                N_low       = N_lower[y + ny * x + ny * nbin * (i - 1)];
                P_low       = P_lower[y + ny * x + ny * nbin * (i - 1)];
                G_pl_low    = G_plus_lower[y + ny * x + ny * nbin * (i - 1)];
                G_min_low   = G_minus_lower[y + ny * x + ny * nbin * (i - 1)];
                g0_low      = g_0_tot_lower[y + ny * x + ny * nbin * (i - 1)];

                // upper part of layer quantities
                w0_up      = w_0_upper[y + ny * x + ny * nbin * (i - 1)];
                del_tau_up = delta_tau_wg_upper[y + ny * x + ny * nbin * (i - 1)];
                M_up       = M_upper[y + ny * x + ny * nbin * (i - 1)];
                N_up       = N_upper[y + ny * x + ny * nbin * (i - 1)];
                P_up       = P_upper[y + ny * x + ny * nbin * (i - 1)];
                G_pl_up    = G_plus_upper[y + ny * x + ny * nbin * (i - 1)];
                G_min_up   = G_minus_upper[y + ny * x + ny * nbin * (i - 1)];
                g0_up      = g_0_tot_upper[y + ny * x + ny * nbin * (i - 1)];

                // improved scattering correction factor E
                E_low = 1.0;
                E_up  = 1.0;

                // improved scattering correction disabled for the following terms -- at least for the moment
                if (scat_corr) {
                    E_up  = E_parameter(w0_up, g0_up, i2s_transition);
                    E_low = E_parameter(w0_low, g0_low, i2s_transition);
                }

                // lower part of layer calculations
                if (del_tau_low < delta_tau_limit) {
                    // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = ((planckband_int[(i - 1) + x * numinterfaces]
                                     + planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                    / 2.0 * (N_low + M_low - P_low));
                }
                else {
                    // non-isothermal solution -- standard case
                    double pgrad_low = (planckband_int[(i - 1) + x * numinterfaces]
                                        - planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                       / del_tau_low;

                    planck_terms =
                        planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * (M_low + N_low)
                        - planckband_int[(i - 1) + x * numinterfaces] * P_low
                        + epsi / (E_low * (1.0 - w0_low * g0_low)) * pgrad_low
                              * (M_low - P_low - N_low);
                }
                flux_terms = P_low * F_up_wg[y + ny * x + ny * nbin * (i - 1)]
                             - N_low * Fc_down_wg[y + ny * x + ny * nbin * (i - 1)];

                direct_terms =
                    Fc_dir_wg[y + ny * x + ny * nbin * (i - 1)] / (-mu_star)
                        * (G_min_low * N_low + G_pl_low * M_low)
                    - F_dir_wg[y + ny * x + ny * nbin * (i - 1)] / (-mu_star) * P_low * G_pl_low;
                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                Fc_up_wg[y + ny * x + ny * nbin * (i - 1)] =
                    1.0 / M_low
                    * (flux_terms
                       + 2.0 * PI * epsi * (1.0 - w0_low) / (E_low - w0_low) * planck_terms
                       + direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (Fc_up_wg[y + ny * x + ny * nbin * (i - 1)] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i - 1,
                               x,
                               y);
                }

                // upper part of layer calculations
                if (del_tau_up < delta_tau_limit) {
                    // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = (planckband_int[i + x * numinterfaces]
                                    + planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                   / 2.0 * (N_up + M_up - P_up);
                }
                else {
                    // non-isothermal solution -- standard case
                    double pgrad_up = (planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)]
                                       - planckband_int[i + x * numinterfaces])
                                      / del_tau_up;

                    planck_terms =
                        planckband_int[i + x * numinterfaces] * (M_up + N_up)
                        - planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * P_up
                        + epsi / (E_up * (1.0 - w0_up * g0_up)) * pgrad_up * (M_up - P_up - N_up);
                }
                flux_terms = P_up * Fc_up_wg[y + ny * x + ny * nbin * (i - 1)]
                             - N_up * F_down_wg[y + ny * x + ny * nbin * i];

                direct_terms =
                    F_dir_wg[y + ny * x + ny * nbin * i] / (-mu_star)
                        * (G_min_up * N_up + G_pl_up * M_up)
                    - Fc_dir_wg[y + ny * x + ny * nbin * (i - 1)] / (-mu_star) * P_up * G_pl_up;

                if (mu_star >= 0.0)
                    direct_terms = 0.0;
                else
                    direct_terms = min(0.0, direct_terms);

                F_up_wg[y + ny * x + ny * nbin * i] =
                    1.0 / M_up
                    * (flux_terms + 2.0 * PI * epsi * (1.0 - w0_up) / (E_up - w0_up) * planck_terms
                       + direct_terms);

                //feedback if flux becomes negative
                if (debug) {
                    if (F_up_wg[y + ny * x + ny * nbin * i] < 0)
                        printf("WARNING WARNING WARNING WARNING -- negative flux found at layer: "
                               "%d, w-index: %d, y-index: %d !!! \n",
                               i,
                               x,
                               y);
                }
            }
        }
    }
}


// calculation of the spectral fluxes, isothermal             case, using thomas algorithm
__global__ void fband_iso_thomas(double* F_down_wg_,      // out
                                 double* F_up_wg_,        // out
                                 double* F_dir_wg_,       // in
                                 double* planckband_lay_, // in
                                 double* w_0_,            // in
                                 double* M_term_,         // in
                                 double* N_term_,         // in
                                 double* P_term_,         // in
                                 double* G_plus_,         // in
                                 double* G_minus_,        // in
                                 double* A_buff_,         // thomas worker
                                 double* B_buff_,         // thomas worker
                                 double* C_buff_,         // thomas worker
                                 double* D_buff_,         // thomas worker
                                 double* C_prime_buff_,   // thomas worker
                                 double* D_prime_buff_,   // thomas worker
                                 double* X_buff_,         // thomas worker
                                 double* g_0_tot_,        // in (clouds)
                                 double* surface_albedo,
                                 bool    singlewalk,
                                 double  Rstar,
                                 double  a,
                                 int     numinterfaces,
                                 int     nbin,
                                 double  f_factor,
                                 double* mu_star_cols,
                                 int     ny,
                                 int     num_cols,
                                 double  epsi,
                                 bool    dir_beam,
                                 bool    clouds,
                                 bool    scat_corr,
                                 bool    debug,
                                 double  i2s_transition) {

    // wavelength bin index
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    // weight index
    int y      = threadIdx.y + blockIdx.y * blockDim.y;
    int c      = blockIdx.z;
    int nlayer = numinterfaces - 1;
    if (x < nbin && y < ny) {
        // Start by applying column offset
        double* F_down_wg      = &(F_down_wg_[c * numinterfaces * nbin * ny]);
        double* F_up_wg        = &(F_up_wg_[c * numinterfaces * nbin * ny]);
        double* F_dir_wg       = &(F_dir_wg_[c * numinterfaces * nbin * ny]);
        double* planckband_lay = &(planckband_lay_[c * (nlayer + 2) * nbin]);
        double* w_0            = &(w_0_[c * nlayer * nbin * ny]);
        double* M_term         = &(M_term_[c * nlayer * nbin * ny]);
        double* N_term         = &(N_term_[c * nlayer * nbin * ny]);
        double* P_term         = &(P_term_[c * nlayer * nbin * ny]);
        double* G_plus         = &(G_plus_[c * nlayer * nbin * ny]);
        double* G_minus        = &(G_minus_[c * nlayer * nbin * ny]);
        double* A_buff         = &(A_buff_[c * numinterfaces * nbin * ny * 4]);
        double* B_buff         = &(B_buff_[c * numinterfaces * nbin * ny * 4]);
        double* C_buff         = &(C_buff_[c * numinterfaces * nbin * ny * 4]);
        double* D_buff         = &(D_buff_[c * numinterfaces * nbin * ny * 2]);
        double* C_prime_buff   = &(C_prime_buff_[c * numinterfaces * nbin * ny * 4]);
        double* D_prime_buff   = &(D_prime_buff_[c * numinterfaces * nbin * ny * 2]);
        double* X_buff         = &(X_buff_[c * numinterfaces * nbin * ny * 2]);
        double* g_0_tot        = &(g_0_tot_[c * nlayer * nbin * ny]);


        // make contiguous address space for worker memory, fastest idx is interface
        // two equations per interface, one matrix block per interface
        int      N       = numinterfaces;
        double4* A       = (double4*)&(A_buff[(x * ny + y) * N * 4]);
        double4* B       = (double4*)&(B_buff[(x * ny + y) * N * 4]);
        double4* C       = (double4*)&(C_buff[(x * ny + y) * N * 4]);
        double2* D       = (double2*)&(D_buff[(x * ny + y) * N * 2]);
        double2* X       = (double2*)&(X_buff[(x * ny + y) * N * 2]);
        double4* C_prime = (double4*)&(C_prime_buff[(x * ny + y) * N * 4]);
        double2* D_prime = (double2*)&(D_prime_buff[(x * ny + y) * N * 2]);

        double mu_star = mu_star_cols[c];

        double F_BOA_up   = 0.0;
        double F_TOA_down = 0.0;

        if (!dir_beam)
            F_TOA_down = f_factor * ((Rstar / a) * (Rstar / a)) * PI
                         * planckband_lay[(N - 1) + x * (numinterfaces - 1 + 2)];

        {

            // BOA boundary -- surface emission

            double w0_N = w_0[y + ny * x + ny * nbin * 0];
            double g0_N = g_0_tot[y + ny * x + ny * nbin * 0];

            // improved scattering correction factor E
            double E_N = 1.0;
            if (scat_corr) {
                E_N = E_parameter(w0_N, g0_N, i2s_transition);
            }

            F_BOA_up = (1 - surface_albedo[x]) * PI * (1.0 - w0_N) / (E_N - w0_N)
                       * planckband_lay[numinterfaces + x * (numinterfaces - 1 + 2)];
            // No upward flux from ghost layer, use stephan boltzman law in net flux computation
            // F_BOA_up = 0.0;
        }

        {
            double w0_0    = w_0[y + ny * x + ny * nbin * 0];
            double M_0     = M_term[y + ny * x + ny * nbin * 0];
            double N_0     = N_term[y + ny * x + ny * nbin * 0];
            double P_0     = P_term[y + ny * x + ny * nbin * 0];
            double G_pl_0  = G_plus[y + ny * x + ny * nbin * 0];
            double G_min_0 = G_minus[y + ny * x + ny * nbin * 0];
            double g0_0    = g_0_tot[y + ny * x + ny * nbin * 0];


            // improved scattering correction factor E
            double E_0 = 1.0;
            if (scat_corr) {
                E_0 = E_parameter(w0_0, g0_0, i2s_transition);
            }

            double B_down = planckband_lay[0 + x * (numinterfaces - 1 + 2)] * (N_0 + M_0 - P_0);

            double I_down =
                min(0.0,
                    -F_dir_wg[y + ny * x + ny * nbin * 0] * (G_min_0 * M_0 + G_pl_0 * N_0)
                        + F_dir_wg[y + ny * x + ny * nbin * (0 + 1)] * P_0 * G_min_0);
            if (mu_star >= 0.0)
                I_down = 0.0;

            // double psi = P;
            // double xi = N;
            // double chi = M;

            A[0].x = 0.0;
            A[0].y = 0.0;
            A[0].z = 0.0;
            A[0].w = 0.0;

            B[0].x = 1.0;
            B[0].y = -1.0 * surface_albedo[x]; //reflected downward beam
            B[0].z = N_0;
            B[0].w = M_0;

            C[0].x = 0.0;
            C[0].y = 0.0;
            C[0].z = 0.0;
            C[0].w = -P_0;

            // BC: upward flux = BOA emission and reflected direct beam
            D[0].x = F_BOA_up + surface_albedo[x] * F_dir_wg[y + ny * x + ny * nbin * 0];
            D[0].y = 2 * PI * epsi * (1.0 - w0_0) / (E_0 - w0_0) * B_down + I_down;
        }

        // fill in inner interfaces
        for (int i = 1; i < N - 1; i++) {
            double w0_down    = w_0[y + ny * x + ny * nbin * i];
            double M_down     = M_term[y + ny * x + ny * nbin * i];
            double N_down     = N_term[y + ny * x + ny * nbin * i];
            double P_down     = P_term[y + ny * x + ny * nbin * i];
            double G_pl_down  = G_plus[y + ny * x + ny * nbin * i];
            double G_min_down = G_minus[y + ny * x + ny * nbin * i];
            double g0_down    = g_0_tot[y + ny * x + ny * nbin * i];

            double w0_up    = w_0[y + ny * x + ny * nbin * (i - 1)];
            double M_up     = M_term[y + ny * x + ny * nbin * (i - 1)];
            double N_up     = N_term[y + ny * x + ny * nbin * (i - 1)];
            double P_up     = P_term[y + ny * x + ny * nbin * (i - 1)];
            double G_pl_up  = G_plus[y + ny * x + ny * nbin * (i - 1)];
            double G_min_up = G_minus[y + ny * x + ny * nbin * (i - 1)];
            double g0_up    = g_0_tot[y + ny * x + ny * nbin * (i - 1)];

            // improved scattering correction factor E
            double E_down = 1.0;
            double E_up   = 1.0;
            if (scat_corr) {
                E_down = E_parameter(w0_down, g0_down, i2s_transition);
                E_up   = E_parameter(w0_up, g0_up, i2s_transition);
            }

            double B_down =
                planckband_lay[i + x * (numinterfaces - 1 + 2)] * (N_down + M_down - P_down);
            double B_up =
                planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * (N_up + M_up - P_up);

            double I_down = min(
                0.0,
                -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min_down * M_down + G_pl_down * N_down)
                    + F_dir_wg[y + ny * x + ny * nbin * (i + 1)] * P_down * G_min_down);

            if (mu_star >= 0.0)
                I_down = 0.0;


            double I_up =
                min(0.0,
                    -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min_up * N_up + G_pl_up * M_up)
                        + F_dir_wg[y + ny * x + ny * nbin * (i - 1)] * P_up * G_pl_up);

            if (mu_star >= 0.0)
                I_up = 0.0;

            // double psi = P;
            // double xi = N;
            // double chi = M;

            A[i].x = -P_up; // -psi
            A[i].y = 0.0;
            A[i].z = 0.0;
            A[i].w = 0.0;

            B[i].x = M_up;   // chi
            B[i].y = N_up;   // xi
            B[i].z = N_down; // xi
            B[i].w = M_down; // chi

            C[i].x = 0.0;
            C[i].y = 0.0;
            C[i].z = 0.0;
            C[i].w = -P_down; // -psi

            D[i].x = 2 * PI * epsi * (1.0 - w0_up) / (E_up - w0_up) * B_up + I_up;
            D[i].y = 2 * PI * epsi * (1.0 - w0_down) / (E_down - w0_down) * B_down + I_down;
        }

        {
            // TOA boundary condition

            double w0_N    = w_0[y + ny * x + ny * nbin * (N - 2)];
            double M_N     = M_term[y + ny * x + ny * nbin * (N - 2)];
            double N_N     = N_term[y + ny * x + ny * nbin * (N - 2)];
            double P_N     = P_term[y + ny * x + ny * nbin * (N - 2)];
            double G_pl_N  = G_plus[y + ny * x + ny * nbin * (N - 2)];
            double G_min_N = G_minus[y + ny * x + ny * nbin * (N - 2)];
            double g0_N    = g_0_tot[y + ny * x + ny * nbin * (N - 2)];

            // improved scattering correction factor E
            double E_N = 1.0;
            if (scat_corr) {
                E_N = E_parameter(w0_N, g0_N, i2s_transition);
            }

            double B_up = planckband_lay[(N - 2) + x * (numinterfaces - 1 + 2)] * (N_N + M_N - P_N);

            double I_up =
                min(0.0,
                    -F_dir_wg[y + ny * x + ny * nbin * (N - 1)] * (G_min_N * N_N + G_pl_N * M_N)
                        + F_dir_wg[y + ny * x + ny * nbin * (N - 2)] * P_N * G_pl_N);

            if (mu_star >= 0.0)
                I_up = 0.0;

            // double psi = P;
            // double xi = N;
            // double chi = M;

            A[N - 1].x = -P_N; // -psi
            A[N - 1].y = 0.0;
            A[N - 1].z = 0.0;
            A[N - 1].w = 0.0;

            B[N - 1].x = M_N; // chi
            B[N - 1].y = N_N; // xi
            B[N - 1].z = 0.0; // xi
            B[N - 1].w = 1.0; // chi

            C[N - 1].x = 0.0;
            C[N - 1].y = 0.0;
            C[N - 1].z = 0.0;
            C[N - 1].w = 0.0; // -psi

            D[N - 1].x = 2 * PI * epsi * (1.0 - w0_N) / (E_N - w0_N) * B_up + I_up;
            D[N - 1].y = F_TOA_down;
        }

        // Solve Thomas algorithm

        thomas_solve(A, B, C, D, C_prime, D_prime, X, N);

        // Fill output data

        for (int i = 0; i < N; i++) {
            // debugging:
            // if (isnan(X[i].x))
            //     printf("F_up[i: % 3d, nbin: % 3d, ny: % 3d] == NaN\n", i, x, y);
            // if (isnan(X[i].y))
            //     printf("F_down[i: % 3d, nbin: % 3d, ny: % 3d] == NaN\n", i, x, y);

            F_up_wg[y + ny * x + ny * nbin * i]   = X[i].x;
            F_down_wg[y + ny * x + ny * nbin * i] = X[i].y;
        }
    }
}
// *****************************************************************************************************


// calculation of the spectral fluxes, non-isothermal case with emphasis on on-the-fly calculations
__global__ void fband_noniso_thomas(double* F_down_wg_,
                                    double* F_up_wg_,
                                    double* F_dir_wg_,
                                    double* Fc_dir_wg_,
                                    double* planckband_lay_,
                                    double* planckband_int_,
                                    double* w_0_upper_,
                                    double* w_0_lower_,
                                    double* delta_tau_wg_upper_,
                                    double* delta_tau_wg_lower_,
                                    double* M_upper_,
                                    double* M_lower_,
                                    double* N_upper_,
                                    double* N_lower_,
                                    double* P_upper_,
                                    double* P_lower_,
                                    double* G_plus_upper_,
                                    double* G_plus_lower_,
                                    double* G_minus_upper_,
                                    double* G_minus_lower_,
                                    double* A_buff_,       // thomas worker
                                    double* B_buff_,       // thomas worker
                                    double* C_buff_,       // thomas worker
                                    double* D_buff_,       // thomas worker
                                    double* C_prime_buff_, // thomas worker
                                    double* D_prime_buff_, // thomas worker
                                    double* X_buff_,       // thomas worker
                                    double* g_0_tot_upper_,
                                    double* g_0_tot_lower_,
                                    double* surface_albedo,
                                    bool    singlewalk,
                                    double  Rstar,
                                    double  a,
                                    int     numinterfaces,
                                    int     nbin,
                                    double  f_factor,
                                    double* mu_star_cols,
                                    int     ny,
                                    int     num_cols,
                                    double  epsi,
                                    double  delta_tau_limit,
                                    bool    dir_beam,
                                    bool    clouds,
                                    bool    scat_corr,
                                    bool    debug,
                                    double  i2s_transition) {
    // wavelength bin index
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    // weight index
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int c = blockIdx.z;

    int nlayer = numinterfaces - 1;
    if (x < nbin && y < ny) {
        double* F_down_wg          = &(F_down_wg_[c * numinterfaces * nbin * ny]);
        double* F_up_wg            = &(F_up_wg_[c * numinterfaces * nbin * ny]);
        double* F_dir_wg           = &(F_dir_wg_[c * numinterfaces * nbin * ny]);
        double* Fc_dir_wg          = &(Fc_dir_wg_[c * numinterfaces * nbin * ny]);
        double* planckband_lay     = &(planckband_lay_[c * (nlayer + 2) * nbin]);
        double* planckband_int     = &(planckband_int_[c * (nlayer + 1) * nbin]);
        double* w_0_upper          = &(w_0_upper_[c * nlayer * nbin * ny]);
        double* w_0_lower          = &(w_0_lower_[c * nlayer * nbin * ny]);
        double* delta_tau_wg_upper = &(delta_tau_wg_upper_[c * nlayer * nbin * ny]);
        double* delta_tau_wg_lower = &(delta_tau_wg_lower_[c * nlayer * nbin * ny]);
        double* M_upper            = &(M_upper_[c * nlayer * nbin * ny]);
        double* M_lower            = &(M_lower_[c * nlayer * nbin * ny]);
        double* N_upper            = &(N_upper_[c * nlayer * nbin * ny]);
        double* N_lower            = &(N_lower_[c * nlayer * nbin * ny]);
        double* P_upper            = &(P_upper_[c * nlayer * nbin * ny]);
        double* P_lower            = &(P_lower_[c * nlayer * nbin * ny]);
        double* G_plus_upper       = &(G_plus_upper_[c * nlayer * nbin * ny]);
        double* G_plus_lower       = &(G_plus_lower_[c * nlayer * nbin * ny]);
        double* G_minus_upper      = &(G_minus_upper_[c * nlayer * nbin * ny]);
        double* G_minus_lower      = &(G_minus_lower_[c * nlayer * nbin * ny]);
        double* A_buff             = &(A_buff_[c * (2 * nlayer + 1) * nbin * ny * 4]);
        double* B_buff             = &(B_buff_[c * (2 * nlayer + 1) * nbin * ny * 4]);
        double* C_buff             = &(C_buff_[c * (2 * nlayer + 1) * nbin * ny * 4]);
        double* D_buff             = &(D_buff_[c * (2 * nlayer + 1) * nbin * ny * 2]);
        double* C_prime_buff       = &(C_prime_buff_[c * (2 * nlayer + 1) * nbin * ny * 4]);
        double* D_prime_buff       = &(D_prime_buff_[c * (2 * nlayer + 1) * nbin * ny * 2]);
        double* X_buff             = &(X_buff_[c * (2 * nlayer + 1) * nbin * ny * 2]);
        double* g_0_tot_upper      = &(g_0_tot_upper_[c * nlayer * nbin * ny]);
        double* g_0_tot_lower      = &(g_0_tot_lower_[c * nlayer * nbin * ny]);

        // make contiguous address space for worker memory, fastest idx is interface
        // two equations per interface, one matrix block per interface
        // Num input layers
        int num_layers = numinterfaces - 1;

        // num interfaces with duplication
        int num_th_layers     = 2 * num_layers;
        int num_th_interfaces = num_th_layers + 1;

        double4* A       = (double4*)&(A_buff[(x * ny + y) * num_th_interfaces * 4]);
        double4* B       = (double4*)&(B_buff[(x * ny + y) * num_th_interfaces * 4]);
        double4* C       = (double4*)&(C_buff[(x * ny + y) * num_th_interfaces * 4]);
        double2* D       = (double2*)&(D_buff[(x * ny + y) * num_th_interfaces * 2]);
        double2* X       = (double2*)&(X_buff[(x * ny + y) * num_th_interfaces * 2]);
        double4* C_prime = (double4*)&(C_prime_buff[(x * ny + y) * num_th_interfaces * 4]);
        double2* D_prime = (double2*)&(D_prime_buff[(x * ny + y) * num_th_interfaces * 2]);

        // Matrix arrays indices and components:
        // 2*i is interfaces
        // 2*i + 1 is layer centers
        // matrix components
        // upward flux: x, y components of matrices (double4)
        // downward flux: z, w component of matrices
        // vectors:
        // upward component: x component of double2
        // downward component: y component of double2

        double w0_up;
        double del_tau_up;
        double M_up;
        double N_up;
        double P_up;
        double G_pl_up;
        double G_min_up;
        double g0_up;
        double E_up;

        double w0_low;
        double del_tau_low;
        double M_low;
        double N_low;
        double P_low;
        double G_pl_low;
        double G_min_low;
        double g0_low;
        double E_low;

        double planck_terms;
        double direct_terms;

        // calculation of downward fluxes factors from TOA to BOA
        for (int i = numinterfaces - 1; i >= 0; i--) {

            // TOA boundary -- incoming stellar flux
            if (i == numinterfaces - 1) {
                double F_TOA = 0.0;
                if (dir_beam)
                    F_TOA = 0.0;
                else
                    F_TOA = f_factor * ((Rstar / a) * (Rstar / a)) * PI
                            * planckband_lay[i + x * (numinterfaces - 1 + 2)];

                // Factors for downward equation, bottom of layer i
                A[2 * (numinterfaces - 1)].z = 0.0;
                A[2 * (numinterfaces - 1)].w = 0.0;

                B[2 * (numinterfaces - 1)].z = 0.0;
                B[2 * (numinterfaces - 1)].w = 1.0;

                C[2 * (numinterfaces - 1)].z = 0.0;
                C[2 * (numinterfaces - 1)].w = 0.0;

                D[2 * (numinterfaces - 1)].y = F_TOA;
            }
            else {
                // upper part of layer quantities
                w0_up      = w_0_upper[y + ny * x + ny * nbin * i];
                del_tau_up = delta_tau_wg_upper[y + ny * x + ny * nbin * i];
                M_up       = M_upper[y + ny * x + ny * nbin * i];
                N_up       = N_upper[y + ny * x + ny * nbin * i];
                P_up       = P_upper[y + ny * x + ny * nbin * i];
                G_pl_up    = G_plus_upper[y + ny * x + ny * nbin * i];
                G_min_up   = G_minus_upper[y + ny * x + ny * nbin * i];
                g0_up      = g_0_tot_upper[y + ny * x + ny * nbin * i];

                // lower part of layer quantities
                w0_low      = w_0_lower[y + ny * x + ny * nbin * i];
                del_tau_low = delta_tau_wg_lower[y + ny * x + ny * nbin * i];
                M_low       = M_lower[y + ny * x + ny * nbin * i];
                N_low       = N_lower[y + ny * x + ny * nbin * i];
                P_low       = P_lower[y + ny * x + ny * nbin * i];
                G_pl_low    = G_plus_lower[y + ny * x + ny * nbin * i];
                G_min_low   = G_minus_lower[y + ny * x + ny * nbin * i];
                g0_low      = g_0_tot_lower[y + ny * x + ny * nbin * i];

                // improved scattering correction factor E
                E_up  = 1.0;
                E_low = 1.0;

                // improved scattering correction disabled for the following terms -- at least for the moment
                if (scat_corr) {
                    E_up  = E_parameter(w0_up, g0_up, i2s_transition);
                    E_low = E_parameter(w0_low, g0_low, i2s_transition);
                }

                // upper part of layer calculations
                if (del_tau_up < delta_tau_limit || P_up == 0.0) {
                    // the isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = (planckband_int[(i + 1) + x * numinterfaces]
                                    + planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                   / 2.0 * (N_up + M_up - P_up);
                }
                else {
                    // the non-isothermal solution -- standard case
                    double pgrad_up = (planckband_lay[i + x * (numinterfaces - 1 + 2)]
                                       - planckband_int[(i + 1) + x * numinterfaces])
                                      / del_tau_up;

                    planck_terms =
                        planckband_lay[i + x * (numinterfaces - 1 + 2)] * (M_up + N_up)
                        - planckband_int[(i + 1) + x * numinterfaces] * P_up
                        + epsi / (E_up * (1.0 - w0_up * g0_up)) * (P_up - M_up + N_up) * pgrad_up;
                }

                direct_terms =
                    -Fc_dir_wg[y + ny * x + ny * nbin * i] * (G_min_up * M_up + G_pl_up * N_up)
                    + F_dir_wg[y + ny * x + ny * nbin * (i + 1)] * G_min_up * P_up;

                direct_terms = min(0.0, direct_terms);

                // Factors for downward equation, center of layer i
                A[2 * i + 1].z = 0.0;
                A[2 * i + 1].w = 0.0;

                B[2 * i + 1].z = N_up;
                B[2 * i + 1].w = M_up;

                C[2 * i + 1].z = 0.0;
                C[2 * i + 1].w = -P_up;

                D[2 * i + 1].y =
                    2.0 * PI * epsi * (1.0 - w0_up) / (E_up - w0_up) * planck_terms + direct_terms;

                // lower part of layer calculations
                if (i == 0) {
                    // Boundary condition for lower layer, equivalent to isotermal solution
                    double pb_lay = planckband_lay[i + x * (numinterfaces - 1 + 2)];
                    planck_terms  = pb_lay * (N_low + M_low - P_low);
                }
                else {
                    {
                        if (del_tau_low < delta_tau_limit || P_low == 0.0) {
                            // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                            planck_terms = (planckband_int[i + x * numinterfaces]
                                            + planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                           / 2.0 * (N_low + M_low - P_low);
                        }
                        else {
                            // non-isothermal solution -- standard case
                            double pgrad_low = (planckband_int[i + x * numinterfaces]
                                                - planckband_lay[i + x * (numinterfaces - 1 + 2)])
                                               / del_tau_low;

                            planck_terms = planckband_int[i + x * numinterfaces] * (M_low + N_low)
                                           - planckband_lay[i + x * (numinterfaces - 1 + 2)] * P_low
                                           + epsi / (E_low * (1.0 - w0_low * g0_low))
                                                 * (P_low - M_low + N_low) * pgrad_low;
                        }
                    }
                }

                direct_terms =
                    -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min_low * M_low + G_pl_low * N_low)
                    + Fc_dir_wg[y + ny * x + ny * nbin * i] * P_low * G_min_low;

                direct_terms = min(0.0, direct_terms);

                // Factors for downward equation, bottom of layer i
                A[2 * i].z = 0.0;
                A[2 * i].w = 0.0;

                B[2 * i].z = N_low;
                B[2 * i].w = M_low;

                C[2 * i].z = 0.0;
                C[2 * i].w = -P_low;

                D[2 * i].y = 2.0 * PI * epsi * (1.0 - w0_low) / (E_low - w0_low) * planck_terms
                             + direct_terms;
            }
        }

        __syncthreads();

        // calculation of upward fluxes from BOA to TOA
        for (int i = 0; i < numinterfaces; i++) {

            // BOA boundary -- surface emission and reflection
            if (i == 0) {
                // No upward flux from ghost layer, use stephan boltzman law in net flux computation
                double F_BOA = 0.0;

                double w0_N = w_0_lower[y + ny * x + ny * nbin * 0];
                double g0_N = g_0_tot_lower[y + ny * x + ny * nbin * 0];

                // improved scattering correction factor E
                double E_N = 1.0;
                if (scat_corr) {
                    E_N = E_parameter(w0_N, g0_N, i2s_transition);
                }

                // printf("%d, %d, %g, %g\n", x, y, w0_N, E_N);
                F_BOA = (1 - surface_albedo[x]) * PI * (1.0 - w0_N) / (E_N - w0_N)
                        * planckband_lay[numinterfaces + x * (numinterfaces - 1 + 2)];

                // Factors for upward equation, top of layer i - bottom of layer i - 1
                A[2 * i].x = 0.0;
                A[2 * i].y = 0.0;

                B[2 * i].x = 1.0;
                B[2 * i].y = -1.0 * surface_albedo[x]; //reflected downward beam

                C[2 * i].x = 0.0;
                C[2 * i].y = 0.0;

                D[2 * i].x = F_BOA + surface_albedo[x] * F_dir_wg[y + ny * x + ny * nbin * i];
            }
            else {
                // lower part of layer quantities
                w0_low      = w_0_lower[y + ny * x + ny * nbin * (i - 1)];
                del_tau_low = delta_tau_wg_lower[y + ny * x + ny * nbin * (i - 1)];
                M_low       = M_lower[y + ny * x + ny * nbin * (i - 1)];
                N_low       = N_lower[y + ny * x + ny * nbin * (i - 1)];
                P_low       = P_lower[y + ny * x + ny * nbin * (i - 1)];
                G_pl_low    = G_plus_lower[y + ny * x + ny * nbin * (i - 1)];
                G_min_low   = G_minus_lower[y + ny * x + ny * nbin * (i - 1)];
                g0_low      = g_0_tot_lower[y + ny * x + ny * nbin * (i - 1)];

                // upper part of layer quantities
                w0_up      = w_0_upper[y + ny * x + ny * nbin * (i - 1)];
                del_tau_up = delta_tau_wg_upper[y + ny * x + ny * nbin * (i - 1)];
                M_up       = M_upper[y + ny * x + ny * nbin * (i - 1)];
                N_up       = N_upper[y + ny * x + ny * nbin * (i - 1)];
                P_up       = P_upper[y + ny * x + ny * nbin * (i - 1)];
                G_pl_up    = G_plus_upper[y + ny * x + ny * nbin * (i - 1)];
                G_min_up   = G_minus_upper[y + ny * x + ny * nbin * (i - 1)];
                g0_up      = g_0_tot_upper[y + ny * x + ny * nbin * (i - 1)];

                // improved scattering correction factor E
                E_low = 1.0;
                E_up  = 1.0;

                // improved scattering correction disabled for the following terms -- at least for the moment
                if (scat_corr) {
                    E_up  = E_parameter(w0_up, g0_up, i2s_transition);
                    E_low = E_parameter(w0_low, g0_low, i2s_transition);
                }

                // lower part of layer calculations
                if (del_tau_low < delta_tau_limit || P_low == 0.0) {
                    // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = ((planckband_int[(i - 1) + x * numinterfaces]
                                     + planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                    / 2.0 * (N_low + M_low - P_low));
                }
                else {
                    // non-isothermal solution -- standard case
                    double pgrad_low = (planckband_int[(i - 1) + x * numinterfaces]
                                        - planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                       / del_tau_low;

                    planck_terms =
                        planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * (M_low + N_low)
                        - planckband_int[(i - 1) + x * numinterfaces] * P_low
                        + epsi / (E_low * (1.0 - w0_low * g0_low)) * pgrad_low
                              * (M_low - P_low - N_low);
                }

                direct_terms = -Fc_dir_wg[y + ny * x + ny * nbin * (i - 1)]
                                   * (G_min_low * N_low + G_pl_low * M_low)
                               + F_dir_wg[y + ny * x + ny * nbin * (i - 1)] * P_low * G_pl_low;

                direct_terms = min(0.0, direct_terms);

                // Factors for upward equation, center of layer i
                A[2 * (i - 1) + 1].x = -P_low;
                A[2 * (i - 1) + 1].y = 0.0;

                B[2 * (i - 1) + 1].x = M_low;
                B[2 * (i - 1) + 1].y = N_low;

                C[2 * (i - 1) + 1].x = 0.0;
                C[2 * (i - 1) + 1].y = 0.0;

                D[2 * (i - 1) + 1].x =
                    2.0 * PI * epsi * (1.0 - w0_low) / (E_low - w0_low) * planck_terms
                    + direct_terms;

                // upper part of layer calculations
                if (del_tau_up < delta_tau_limit || P_up == 0.0) {
                    // isothermal solution -- taken if optical depth so small that numerical instabilities may occur
                    planck_terms = (planckband_int[i + x * numinterfaces]
                                    + planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)])
                                   / 2.0 * (N_up + M_up - P_up);
                }
                else {
                    // non-isothermal solution -- standard case
                    double pgrad_up = (planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)]
                                       - planckband_int[i + x * numinterfaces])
                                      / del_tau_up;

                    planck_terms =
                        planckband_int[i + x * numinterfaces] * (M_up + N_up)
                        - planckband_lay[(i - 1) + x * (numinterfaces - 1 + 2)] * P_up
                        + epsi / (E_up * (1.0 - w0_up * g0_up)) * pgrad_up * (M_up - P_up - N_up);
                }

                direct_terms =
                    -F_dir_wg[y + ny * x + ny * nbin * i] * (G_min_up * N_up + G_pl_up * M_up)
                    + Fc_dir_wg[y + ny * x + ny * nbin * (i - 1)] * P_up * G_pl_up;

                direct_terms = min(0.0, direct_terms);

                // Factors for upward equation, top of layer i - bottom of layer i - 1
                A[2 * i].x = -P_up;
                A[2 * i].y = 0.0;

                B[2 * i].x = M_up;
                B[2 * i].y = N_up;

                C[2 * i].x = 0.0;
                C[2 * i].y = 0.0;

                D[2 * i].x =
                    2.0 * PI * epsi * (1.0 - w0_up) / (E_up - w0_up) * planck_terms + direct_terms;
            }
        }
        // Solve Thomas algorithm

        thomas_solve(A, B, C, D, C_prime, D_prime, X, num_th_interfaces);

        // Fill output data

        for (int i = 0; i < numinterfaces; i++) {
            // if (isnan(X[i].x))
            //     printf("F_up[i: % 3d, nbin: % 3d, ny: % 3d] == NaN\n", i, x, y);
            // if (isnan(X[i].y))
            //     printf("F_down[i: % 3d, nbin: % 3d, ny: % 3d] == NaN\n", i, x, y);

            F_up_wg[y + ny * x + ny * nbin * i]   = X[2 * i].x;
            F_down_wg[y + ny * x + ny * nbin * i] = X[2 * i].y;
        }
    }
}
