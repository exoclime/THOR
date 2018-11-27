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
// Description: Deep hot jupiter benchmark test (Heng+ 2011)
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
//       [2] Heng et al. 2011
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

__global__ void deepHJ(double *Mh_d,
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
        //double sigma;
        double pre;
        //double psm1;
        double lat = lonlat_d[id * 2 + 1];
        double lon = lonlat_d[id * 2];
        if (lon < 0) lon += 2 * M_PI;

        // Deeeeeep Hot Jupiter stuff....

        double Teq_hs;
        double kv_hs, kt_hs;

        double Tnight, Tday;
        double Ptil;
        double logtrad;

        ////////////
        //      Calculates surface pressure
        //psm1 = pressure_d[id*nv + 1] - Rho_d[id*nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
        //ps = 0.5*(pressure_d[id*nv + 0] + psm1);

        pre = pressure_d[id * nv + lev];

        //sigma = (pre/ps);
        Ptil = log10(pre / 100000); // log of pressure in bars

        if (pre <= 1e6) { //pressure less than ten bars
            Tnight = 1388.2145 + 267.66586 * Ptil - 215.53357 * Ptil * Ptil
                     + 61.814807 * Ptil * Ptil * Ptil + 135.68661 * Ptil * Ptil * Ptil * Ptil
                     + 2.01149044 * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 40.907246 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 19.015628 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 3.8771634 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 0.38413901 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 0.015089084 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
            Tday = 2149.9581 + 4.1395571 * Ptil - 186.24851 * Ptil * Ptil
                   + 135.52524 * Ptil * Ptil * Ptil + 106.20433 * Ptil * Ptil * Ptil * Ptil
                   - 35.851966 * Ptil * Ptil * Ptil * Ptil * Ptil
                   - 50.022826 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                   - 18.462489 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                   - 3.3319965 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                   - 0.30295925 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                   - 0.011122316 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
            logtrad = 5.4659686 + 1.4940124 * Ptil + 0.66079196 * Ptil * Ptil
                      + 0.16475329 * Ptil * Ptil * Ptil + 0.014241552 * Ptil * Ptil * Ptil * Ptil;
            kt_hs = pow(10, -logtrad); //      Temperature forcing constant.
        }
        else {
            /*Tnight = 5529.7168 - 6869.6504*Ptil + 4142.7231*Ptil*Ptil \
                   - 936.23053*Ptil*Ptil*Ptil + 87.120975*Ptil*Ptil*Ptil*Ptil;*/
            //The above seems to be a very poor fit to Tiro... switch to eq B4 instead
            Tnight = 1696.6986 + 132.2318 * Ptil - 174.30459 * Ptil * Ptil
                     + 12.579612 * Ptil * Ptil * Ptil + 59.513639 * Ptil * Ptil * Ptil * Ptil
                     + 9.6706522 * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 4.1136048 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     - 1.0632301 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     + 0.064400203 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     + 0.035974396 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
                     + 0.0025740066 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
            Tday  = Tnight;
            kt_hs = 0.0; //      Temperature forcing constant.
        }

        // Calculate Teq from eqn 26 (Heng, Menou, and Phillips 2011)
        if ((lon >= M_PI / 2) && (lon < 3 * M_PI / 2)) {
            Teq_hs = pow(Tnight * Tnight * Tnight * Tnight + (Tday * Tday * Tday * Tday - Tnight * Tnight * Tnight * Tnight) * cos(lon - M_PI) * cos(lat), 0.25);
        }
        else {
            Teq_hs = Tnight;
        }

        //      Momentum dissipation constant.
        kv_hs = 0.0; //no boundary layer friction

        //      Update momenta
        for (int k = 0; k < 3; k++) Mh_d[id * 3 * nv + lev * 3 + k] = Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + kv_hs * time_step);
        ;
        //      Update temperature
        temperature_d[id * nv + lev] -= kt_hs * time_step * (temperature_d[id * nv + lev] - Teq_hs);
    }
}
