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
// Description: Earth benchmark test.
//
//
// Method: The temperature is forced using a Newtonian cooling code, and the
//         boundary layer is represented by a Rayleigh friction scheme.
//
// Known limitations: None.
//
// Known issues: None.
//
// If you use this code please cite the following references:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
//       [2] Held, I. M., & Suarez, M. J. 1994, Bullentin of the American
//           Meteorological Society
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

__global__ void held_suarez(double *Mh_d,
                            double *pressure_d,
                            double *Rho_d,
                            double *temperature_d,
                            double  Gravit,
                            double  Cp,
                            double  Rd,
                            double *Altitude_d,
                            double *Altitudeh_d,
                            double *lonlat_d,
                            double  time_step,
                            int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {

        //      Parameters for the forcing and dissipation
        double sigma;
        double sigma0;
        double sigmab  = 0.7;
        double ka      = (1.0 / 40.0) * (1.0 / 86400.0);
        double kf      = 1.0 / 86400.0;
        double ks      = (1.0 / 4.0) * (1.0 / 86400.0);
        double kappa   = Rd / Cp;
        double dTy     = 60.0;
        double dthetaz = 10.0;
        double ps, pre;
        double psm1;
        double p0  = 100000.0;
        double lat = lonlat_d[id * 2 + 1];

        double Teq_hs;
        double kv_hs, kt_hs;
        ////////////
        //      Calculates surface pressure
        psm1 = pressure_d[id * nv + 1]
               - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
        ps = 0.5 * (pressure_d[id * nv + 0] + psm1);

        pre = pressure_d[id * nv + lev];

        sigma  = (pre / ps);
        sigma0 = (pre / p0);

        //      Equilibrium temperature.
        Teq_hs = max(200.0,
                     (315.0 - dTy * pow(sin(lat), 2.0) - dthetaz * log(sigma0) * pow(cos(lat), 2.0))
                         * pow(sigma0, kappa));

        //      Temperature forcing constant.
        kt_hs = ka + (ks - ka) * max(0.0, (sigma - sigmab) / (1.0 - sigmab)) * pow(cos(lat), 4.0);

        //      Momentum dissipation constant.
        kv_hs = kf * max(0.0, (sigma - sigmab) / (1.0 - sigmab));

        //      Update momenta
        for (int k = 0; k < 3; k++)
            Mh_d[id * 3 * nv + lev * 3 + k] =
                Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + kv_hs * time_step);

        //      Update temperature
        temperature_d[id * nv + lev] -= kt_hs * time_step * (temperature_d[id * nv + lev] - Teq_hs);
    }
}
