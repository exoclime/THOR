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
// Description: Shallow hot jupiter benchmark test.
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
//       [2] Heng et al. 2011, Rauscher & Menou ?
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

__global__ void shallowHJ(double *Mh_d,
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
        //double sigmab = 0.7;
        //double ka     = (1.0/40.0) * (1.0/86400.0);
        //double ks     = (1.0/4.0) * (1.0/86400.0);
        double dTy = 200.0;
        double ps, pre;
        double psm1;
        double lat = lonlat_d[id * 2 + 1];
        double lon = lonlat_d[id * 2];

        // Hot Jupiter stuff....
        double zstra       = 2e6;
        double sigma_stra  = 0.12;
        double Gamma_trop  = 2e-4;
        double T_surf      = 1600;
        double Delta_Tstra = 10;
        double beta_trop;
        double T_vert;

        double Teq_hs;
        double kv_hs, kt_hs;
        ////////////
        //      Calculates surface pressure
        psm1 = pressure_d[id * nv + 1]
               - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
        ps = 0.5 * (pressure_d[id * nv + 0] + psm1);

        pre = pressure_d[id * nv + lev];

        sigma = (pre / ps);

        // Calculate T_vert from Equation 23 & 24 (Heng, Menou, and Phillips 2011)
        if (Altitude_d[lev] <= zstra) {
            beta_trop = sin(0.5 * M_PI * (sigma - sigma_stra) / (1.0 - sigma_stra));
            T_vert    = T_surf - Gamma_trop * (zstra + 0.5 * (Altitude_d[lev] - zstra))
                     + sqrt(pow(0.5 * Gamma_trop * (Altitude_d[lev] - zstra), 2)
                            + Delta_Tstra * Delta_Tstra);
        }
        else {
            beta_trop = 0.0;
            T_vert    = T_surf - Gamma_trop * zstra + Delta_Tstra;
        }

        //      Equilibrium temperature. Equation 25 from Heng, Menou, and Phillips 2011
        Teq_hs = T_vert + beta_trop * dTy * cos(lon - M_PI) * cos(lat);

        //      Temperature forcing constant.
        kt_hs = 1 / (1.5e5);

        //      Momentum dissipation constant.
        kv_hs = 0.0; //no boundary layer friction

        //      Update momentum
        for (int k = 0; k < 3; k++)
            Mh_d[id * 3 * nv + lev * 3 + k] =
                Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + kv_hs * time_step);
        ;
        //      Update temperature
        temperature_d[id * nv + lev] -= kt_hs * time_step * (temperature_d[id * nv + lev] - Teq_hs);
    }
}
