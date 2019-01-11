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
// Description: Test of acoustic wave behavior
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

__global__ void acoustic_test(double *pressure_d,
                              double *Rho_d,
                              double *temperature_d,
                              double  Rd,
                              double *Altitude_d,
                              double *lonlat_d,
                              double  Top_altitude,
                              int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double dp, R, lambda0, phi0, vmode, r, g, f;
        double lat = lonlat_d[id * 2 + 1];
        double lon = lonlat_d[id * 2];

        vmode   = 1;         // vertical mode
        dp      = 100;       // pressure perturbation (Pa)
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

        pressure_d[id * nv + lev] += dp * f * g;
        Rho_d[id * nv + lev] = pressure_d[id * nv + lev] / Rd / temperature_d[id * nv + lev];
    }
}
