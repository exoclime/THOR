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

#include <math.h>

__global__ void dry_conv_adj(double *Pressure_d,    // Pressure [Pa]
                             double *Pressureh_d,   // Mid-point pressure [Pa]
                             double *Temperature_d, // Temperature [K]
                             double *pt_d,          // Potential temperature [K]
                             double *Rho_d,         // Density [m^3/kg]
                             double  Cp,            // Specific heat capacity [J/kg/K]
                             double  Rd,            // Gas constant [J/kg/K]
                             double  Gravit,        // Gravity [m/s^2]
                             double *Altitude_d,    // Altitudes of the layers
                             double *Altitudeh_d,   // Altitudes of the interfaces
                             int     num,           // Number of columns
                             int     nv) {              // Vertical levels
    //
    //  Description: Mixes entropy vertically on statically unstable columns
    //

    int         id = blockIdx.x * blockDim.x + threadIdx.x;

    // Local arrays
    // double      theta[nv];
    // double      ph[nv + 1];

    // stability threshold
    double      stable = 0.0;

    double      ps, psm;
    double      pp, ptop;

    double      xi, xip, xim, a, b;

    if (id < num) {
        // calculate pressure at the interfaces
        for (int lev = 0; lev <= nv; lev++) {
            if (lev == 0) {
                // extrapolate to lower boundary
                psm = Pressure_d[id * nv + 1]
                      - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                ps                             = 0.5 * (Pressure_d[id * nv + 0] + psm);
                Pressureh_d[id * (nv + 1) + 0] = ps;
            }
            else if (lev == nv) {
                // extrapolate to top boundary
                pp = Pressure_d[id * nv + nv - 2]
                     - Rho_d[id * nv + nv - 1] * Gravit
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0; //prevents pressure from going negative
                ptop                             = 0.5 * (Pressure_d[id * nv + nv - 1] + pp);
                Pressureh_d[id * (nv + 1) + lev] = ptop;
            }
            else {
                // interpolation between layers
                xi  = Altitudeh_d[lev];
                xim = Altitude_d[lev - 1];
                xip = Altitude_d[lev];
                a   = (xi - xip) / (xim - xip);
                b   = (xi - xim) / (xip - xim);
                Pressureh_d[id * (nv + 1) + lev] =
                    Pressure_d[id * nv + lev - 1] * a + Pressure_d[id * nv + lev] * b;
            }
        }

        // Compute Potential Temperature
        for (int lev = 0; lev < nv; lev++) {
            pt_d[id * nv + lev] =
                Temperature_d[id * nv + lev]
                * pow(Pressureh_d[id * (nv + 1) + 0] / Pressure_d[id * nv + lev], Rd / Cp);
        }

        bool done_col = false;
        while (done_col == false) { // Unstable  column?
            int top = 0;
            int bot = nv - 1;

            for (int lev = 0; lev < nv - 1; lev++) {
                // sweep upward, find lowest unstable layer
                if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] < stable) {
                    if (bot > lev)
                        bot = lev;
                }
            }

            for (int lev = bot; lev < nv - 1; lev++) {
                // sweep upward from unstable layer, find top
                if (pt_d[id * nv + lev + 1] - pt_d[id * nv + lev] > stable) {
                    top = lev;
                    break;
                }
                else {
                    top = nv - 1;
                }
            }

            if (bot < nv - 1) {
                int    extend = 1;
                double thnew;

                while (extend == 1) {
                    double h   = 0.0; //Enthalpy;
                    double sum = 0.0;
                    extend     = 0;

                    for (int lev = bot; lev <= top; lev++) {
                        // calc adiabatic pressure, integrate upward for new pot. temp.
                        double pu = Pressureh_d[id * (nv + 1) + lev + 1];
                        double pl = Pressureh_d[id * (nv + 1) + lev];
                        double pi = pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0],
                                        Rd / Cp); // adiabatic pressure wrt bottom of column
                        double deltap = pl - pu;

                        h   = h + pt_d[id * nv + lev] * pi * deltap;
                        sum = sum + pi * deltap;
                    }
                    thnew = h / sum;

                    // if (bot <= 0 && top >= nv - 1) {
                    //     // no need to extend again
                    //     extend = 0;
                    // }

                    if (bot > 0) {
                        // repeat if new pot. temp. is less than lower boundary p.t.
                        if ((thnew - pt_d[id * nv + bot - 1]) < stable) {
                            bot    = bot - 1;
                            extend = 1;
                        }
                    }

                    if (top < nv - 1) {
                        // repeat if new pot. temp. is greater p.t. above
                        if ((pt_d[id * nv + top + 1] - thnew) < stable) {
                            top    = top + 1;
                            extend = 1;
                        }
                    }
                }

                for (int lev = bot; lev <= top; lev++) {
                    pt_d[id * nv + lev] = thnew; // set new potential temperature
                }
            }
            else {
                done_col = true; //no unstable layers
            }
        }

        // Compute Temperature & pressure from potential temperature
        for (int lev = 0; lev < nv; lev++) {
            Temperature_d[id * nv + lev] =
                pt_d[id * nv + lev]
                * pow(Pressure_d[id * nv + lev] / Pressureh_d[id * (nv + 1) + 0], Rd / Cp);
            Pressure_d[id * nv + lev] = Temperature_d[id * nv + lev] * Rd * Rho_d[id * nv + lev];
        }
    }
}
