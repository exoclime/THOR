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
// Description: Test of gravity wave behavior
//
//
// Method:
//
// Known limitations: None.
//
// Known issues: None.
//
// If you use this code please cite the following references:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
//       [2] Tomita & Satoh 2004
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

void setup_gwave_TP(double *pressure_h,
                    double *temperature_h,
                    double *Altitude_h,
                    int     nv,
                    int     id,
                    int     lev) {

    // for now, assume N = 0.01 as in case 1 of T&S
    // RD: i don't think i'm using this function (should prob remove it)
    double logP, T;
    double z = Altitude_h[lev];

    // polyfits for p and T
    logP = -4.30649701e-43 * pow(z, 10) + 6.96539355e-38 * pow(z, 9) - 4.82105798e-33 * pow(z, 8)
           + 1.86013495e-28 * pow(z, 7) - 4.37180502e-24 * pow(z, 6) + 6.42032008e-20 * pow(z, 5)
           - 5.81443216e-16 * pow(z, 4) + 3.05991707e-12 * pow(z, 3) - 8.91875000e-09 * pow(z, 2)
           - 4.07764528e-05 * z + 4.99994103;

    pressure_h[id * nv + lev] = pow(10, logP);

    T = 4.24534875e-43 * pow(z, 10) - 8.20325495e-38 * pow(z, 9) + 6.66627493e-33 * pow(z, 8)
        - 2.97343975e-28 * pow(z, 7) + 7.96256373e-24 * pow(z, 6) - 1.31454761e-19 * pow(z, 5)
        + 1.32002272e-15 * pow(z, 4) - 7.75487103e-12 * pow(z, 3) - 1.15941319e-08 * pow(z, 2)
        - 6.71406985e-03 * z + 3.00000267e+02;
    temperature_h[id * nv + lev] = T;
}

__global__ void gwave_test(double *pressure_d,
                           double *Rho_d,
                           double *temperature_d,
                           double  Rd,
                           double  Cp,
                           double  P_Ref,
                           double *Altitude_d,
                           double *lonlat_d,
                           double  Top_altitude,
                           int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double dpt, R, lambda0, phi0, vmode, r, g, f;
        double lat   = lonlat_d[id * 2 + 1];
        double lon   = lonlat_d[id * 2];
        double kappa = Rd / Cp, pt;

        vmode   = 1;         // vertical mode
        dpt     = 10;        // potential temp perturbation (K)
        R       = 1.0 / 3.0; // distance cutoff of perturbation
        lambda0 = 0;         //longitude of perturbation
        phi0    = 0;         //latitude of perturbation
        r       = acos(sin(phi0) * sin(lat) + cos(phi0) * cos(lat) * cos(lon - lambda0));
        g       = sin(vmode * M_PI * Altitude_d[lev] / Top_altitude);
        if (r < R) {
            f = 0.5 * (1 + cos(M_PI * r / R));
        }
        else {
            f = 0.0;
        }

        pt = temperature_d[id * nv + lev] * pow(pressure_d[id * nv + lev] / P_Ref, -kappa);

        pt += dpt * f * g; // apply perturbation to potential temperature
        temperature_d[id * nv + lev] = pt * pow(pressure_d[id * nv + lev] / P_Ref, kappa);
        Rho_d[id * nv + lev] = pressure_d[id * nv + lev] / (Rd * temperature_d[id * nv + lev]);
    }
}
