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
// Kernels for value interpolation from lookup tables
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

#include "interpolate_values.h"

#include "physics_constants.h"

// temperature interpolation for the non-isothermal layers
__global__ void interpolate_temperature(double* tlay, double* tint, int numinterfaces) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    int c = blockIdx.z;

    if (0 < i && i < numinterfaces - 1) {
        tint[i + c * numinterfaces] =
            tlay[i - 1 + c * numinterfaces]
            + 0.5 * (tlay[i + c * numinterfaces] - tlay[i - 1 + c * numinterfaces]);
    }
    if (i == 0) {
        tint[i + c * numinterfaces] =
            tlay[numinterfaces - 1 + c * numinterfaces]; // set equal to the surface/BOA temperature
    }
    if (i == numinterfaces - 1) {
        tint[i + c * numinterfaces] =
            tlay[i - 1 + c * numinterfaces]
            + 0.5 * (tlay[i - 1 + c * numinterfaces] - tlay[i - 2 + c * numinterfaces]);
    }
}

// interpolates the Planck function for the layer temperatures from the pre-tabulated values
__global__ void planck_interpol_layer(double* temp,           // in
                                      double* planckband_lay, // out
                                      double* planck_grid,    // in
                                      double* starflux,       // in
                                      bool    null_planck_function,
                                      bool    realstar,
                                      int     numlayers,
                                      int     nwave,
                                      int     dim,
                                      int     step) {

    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int c = blockIdx.z;

    if (x < nwave && i < numlayers + 2) {

        planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] = 0.0;

        // getting the stellar flux --- is redundant to do it every interpolation, but probably has negligible costs ...
        if (i == numlayers) {
            if (realstar) {
                planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] =
                    starflux[x] / PI;
            }
            else {
                planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] =
                    planck_grid[x + dim * nwave];
            }
        }
        else {
            double t;

            // interpolating for layer temperatures
            if (i < numlayers) {
                t = (temp[i + c * (numlayers + 1)] - 1.0) / step;
            }
            // interpolating for below (surface/BOA) temperature
            if (i == numlayers + 1) {
                t = (temp[numlayers + c * (numlayers + 1)] - 1.0) / step;
                t = t + 1 - 1;
            }

            t = max(0.001, min(dim - 1.001, t));

            int tdown = floor(t);
            int tup   = ceil(t);

            if (tdown != tup) {
                planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] =
                    planck_grid[x + tdown * nwave] * (tup - t)
                    + planck_grid[x + tup * nwave] * (t - tdown);
            }
            if (tdown == tup) {
                planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] =
                    planck_grid[x + tdown * nwave];
            }
            if (null_planck_function)
                planckband_lay[i + x * (numlayers + 2) + c * (numlayers + 2) * nwave] = 0.0;
        }
    }
}

// interpolates the Planck function for the interface temperatures from the pre-tabulated values
__global__ void planck_interpol_interface(double* temp,           // in
                                          double* planckband_int, // out
                                          double* planck_grid,    // in
                                          bool    null_planck_function,
                                          int     numinterfaces,
                                          int     nwave,
                                          int     dim,
                                          int     step) {

    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int c = blockIdx.z;

    if (x < nwave && i < numinterfaces) {

        double t = (temp[i + c * numinterfaces] - 1.0) / step;

        t = max(0.001, min(dim - 1.001, t));

        int tdown = floor(t);
        int tup   = ceil(t);

        if (tdown != tup) {
            planckband_int[i + x * numinterfaces + c * numinterfaces * nwave] =
                planck_grid[x + tdown * nwave] * (tup - t)
                + planck_grid[x + tup * nwave] * (t - tdown);
        }
        if (tdown == tup) {
            planckband_int[i + x * numinterfaces + c * numinterfaces * nwave] =
                planck_grid[x + tdown * nwave];
        }
        if (null_planck_function)
            planckband_int[i + x * numinterfaces + c * numinterfaces * nwave] = 0.0;
    }
}


// interpolate layer and interface opacities from opacity table
__global__ void interpolate_opacities(
    double* temp,       // in, layer temperature
    double* opactemp,   // in, opac reference table temperatures
    double* press,      // in, layer pressure
    double* opacpress,  // in, opac reference table pressure
    double* ktable,     // in, reference opacity table
    double* opac,       // out, opacities
    double* crosstable, // in, reference scatter cross-section
    double* scat_cross, // out, scattering cross-sections
    int     npress,     // in, ref table pressures count
    int     ntemp,      // in, ref table temperatures count
    int     ny,         // in, ref table gaussians y-points count
    int     nbin,       // in, ref wavelength bin count
    double  opaclimit,  // in, opacity limit for max cutoff for low wavelength bin idx
    int     temp_num_per_col,
    int     press_num_per_col,
    int     nlay_or_nint // in, number of layer or interfaces (physical position bins count)
) {

    int x = threadIdx.x + blockIdx.x * blockDim.x; // wavelength bin index
    int i = threadIdx.y + blockIdx.y * blockDim.y; // volume element (layer or interface) bin index
    int c = blockIdx.z;

    if (x < nbin && i < nlay_or_nint) {

        // looks like it's used to clip the bottom of interpolated opacity to opac_limit when below a certain
        // wavelength threshold
        // see comment above in parameters. "opacity limit for max cutoff for low wavelength bin idx"
        int x_1micron = lrint(nbin * 2.0 / 3.0);

        double deltaopactemp = (opactemp[ntemp - 1] - opactemp[0]) / (ntemp - 1.0);
        double deltaopacpress =
            (log10(opacpress[npress - 1]) - log10(opacpress[0])) / (npress - 1.0);
        double t = (temp[i + c * (temp_num_per_col)] - opactemp[0]) / deltaopactemp;

        t = min(ntemp - 1.001, max(0.001, t));

        int tdown = floor(t);
        int tup   = ceil(t);

        double p = (log10(press[i + c * press_num_per_col]) - log10(opacpress[0])) / deltaopacpress;


        // do the cloud deck
        // This comes from helios, not used in alfrodull, clouds are handled later
        double k_cloud = 0.0; //1e-1 * norm_pdf(log10(press[i]),0,1);

        p = min(npress - 1.001, max(0.001, p));

        int pdown = floor(p);
        int pup   = ceil(p);

        if (pdown != pup && tdown != tup) {
            for (int y = 0; y < ny; y++) {
                double interpolated_opac =
                    ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tdown] * (pup - p)
                        * (tup - t)
                    + ktable[y + ny * x + ny * nbin * pup + ny * nbin * npress * tdown]
                          * (p - pdown) * (tup - t)
                    + ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tup] * (pup - p)
                          * (t - tdown)
                    + ktable[y + ny * x + ny * nbin * pup + ny * nbin * npress * tup] * (p - pdown)
                          * (t - tdown);

                if (x < x_1micron) {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        max(interpolated_opac, opaclimit);
                }
                else {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        interpolated_opac;
                }

                opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] += k_cloud;
            }

            scat_cross[x + nbin * i + c * nlay_or_nint * nbin] =
                crosstable[x + nbin * pdown + nbin * npress * tdown] * (pup - p) * (tup - t)
                + crosstable[x + nbin * pup + nbin * npress * tdown] * (p - pdown) * (tup - t)
                + crosstable[x + nbin * pdown + nbin * npress * tup] * (pup - p) * (t - tdown)
                + crosstable[x + nbin * pup + nbin * npress * tup] * (p - pdown) * (t - tdown);
        }

        if (tdown == tup && pdown != pup) {
            for (int y = 0; y < ny; y++) {
                double interpolated_opac =
                    ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tdown] * (pup - p)
                    + ktable[y + ny * x + ny * nbin * pup + ny * nbin * npress * tdown]
                          * (p - pdown);
                if (x < x_1micron) {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        max(interpolated_opac, opaclimit);
                }
                else {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        interpolated_opac;
                }

                opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] += k_cloud;
            }

            scat_cross[x + nbin * i + c * nlay_or_nint * nbin] =
                crosstable[x + nbin * pdown + nbin * npress * tdown] * (pup - p)
                + crosstable[x + nbin * pup + nbin * npress * tdown] * (p - pdown);
        }

        if (pdown == pup && tdown != tup) {
            for (int y = 0; y < ny; y++) {
                double interpolated_opac =
                    ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tdown] * (tup - t)
                    + ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tup]
                          * (t - tdown);
                if (x < x_1micron) {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        max(interpolated_opac, opaclimit);
                }
                else {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        interpolated_opac;
                }

                opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] += k_cloud;
            }

            scat_cross[x + nbin * i + c * nlay_or_nint * nbin] =
                crosstable[x + nbin * pdown + nbin * npress * tdown] * (tup - t)
                + crosstable[x + nbin * pdown + nbin * npress * tup] * (t - tdown);
        }

        if (tdown == tup && pdown == pup) {
            for (int y = 0; y < ny; y++) {

                double interpolated_opac =
                    ktable[y + ny * x + ny * nbin * pdown + ny * nbin * npress * tdown];

                if (x < x_1micron) {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        max(interpolated_opac, opaclimit);
                }
                else {
                    opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] =
                        interpolated_opac;
                }

                opac[y + ny * x + ny * nbin * i + c * nlay_or_nint * ny * nbin] += k_cloud;
            }

            scat_cross[x + nbin * i + c * nlay_or_nint * nbin] =
                crosstable[x + nbin * pdown + nbin * npress * tdown];
        }
    }
}


// interpolate the mean molecular mass for each layer
__global__ void meanmolmass_interpol(double* temp,          // in
                                     double* opactemp,      // in
                                     double* meanmolmass,   // out
                                     double* opac_meanmass, // in
                                     double* press,         // in
                                     double* opacpress,     // in
                                     int     npress,
                                     int     ntemp,
                                     int     temp_num_per_col,
                                     int     press_num_per_col,
                                     int     ninterface) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int c = blockIdx.z;

    if (i < ninterface) {

        double deltaopactemp = (opactemp[ntemp - 1] - opactemp[0]) / (ntemp - 1.0);
        double deltaopacpress =
            (log10(opacpress[npress - 1]) - log10(opacpress[0])) / (npress - 1.0);
        double t = (temp[i + c * temp_num_per_col] - opactemp[0]) / deltaopactemp;

        t = min(ntemp - 1.001, max(0.001, t));

        int tdown = floor(t);
        int tup   = ceil(t);

        double p = (log10(press[i + c * press_num_per_col]) - log10(opacpress[0])) / deltaopacpress;

        p = min(npress - 1.001, max(0.001, p));

        int pdown = floor(p);
        int pup   = ceil(p);

        if (tdown != tup && pdown != pup) {
            meanmolmass[i + c * ninterface] =
                opac_meanmass[pdown + npress * tdown] * (pup - p) * (tup - t)
                + opac_meanmass[pup + npress * tdown] * (p - pdown) * (tup - t)
                + opac_meanmass[pdown + npress * tup] * (pup - p) * (t - tdown)
                + opac_meanmass[pup + npress * tup] * (p - pdown) * (t - tdown);
        }
        if (tdown != tup && pdown == pup) {
            meanmolmass[i + c * ninterface] = opac_meanmass[pdown + npress * tdown] * (tup - t)
                                              + opac_meanmass[pdown + npress * tup] * (t - tdown);
        }
        if (tdown == tup && pdown != pup) {
            meanmolmass[i + c * ninterface] = opac_meanmass[pdown + npress * tdown] * (pup - p)
                                              + opac_meanmass[pup + npress * tdown] * (p - pdown);
        }
        if (tdown == tup && pdown == pup) {
            meanmolmass[i + c * ninterface] = opac_meanmass[pdown + npress * tdown];
        }
    }
}
