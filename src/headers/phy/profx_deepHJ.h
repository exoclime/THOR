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
//       [2] Held, I. M., & Suarez, M. J. 1994, Bullentin of the American
//           Meteorological Society
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

__global__ void deepHJ(double *Mh_d,
                       double *pressure_d,
                       double *Rho_d,
                       double *temperature_d,
                       double *profx_dP_d,
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
        if (lon < 0)
            lon += 2 * M_PI;

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
        if (pre >= 100) {
            Ptil = log10(pre / 100000); // log of pressure in bars
        }
        else { // here we are outside the bounds of the poly fits, make isothermal
            Ptil = log10(100.0 / 100000);
        }


        if (pre <= 1e6) { //pressure less than ten bars
            // Tnight =
            //     1388.2145 + 267.66586 * Ptil - 215.53357 * Ptil * Ptil
            //     + 61.814807 * Ptil * Ptil * Ptil + 135.68661 * Ptil * Ptil * Ptil * Ptil
            //     + 2.01149044 * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 40.907246 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 19.015628 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 3.8771634 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 0.38413901 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 0.015089084 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
            // Tday =
            //     2149.9581 + 4.1395571 * Ptil - 186.24851 * Ptil * Ptil
            //     + 135.52524 * Ptil * Ptil * Ptil + 106.20433 * Ptil * Ptil * Ptil * Ptil
            //     - 35.851966 * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 50.022826 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 18.462489 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 3.3319965 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 0.30295925 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //     - 0.011122316 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil;
            logtrad = 5.4659686 + 1.4940124 * Ptil + 0.66079196 * Ptil * Ptil
                      + 0.16475329 * Ptil * Ptil * Ptil + 0.014241552 * Ptil * Ptil * Ptil * Ptil;
            kt_hs = pow(10, -logtrad); //      Temperature forcing constant.
        }
        else {
            /*Tnight = 5529.7168 - 6869.6504*Ptil + 4142.7231*Ptil*Ptil \
                   - 936.23053*Ptil*Ptil*Ptil + 87.120975*Ptil*Ptil*Ptil*Ptil;*/
            //The above seems to be a very poor fit to Tiro... switch to eq B4 instead
            // Tnight = 1696.6986 + 132.2318 * Ptil - 174.30459 * Ptil * Ptil
            //          + 12.579612 * Ptil * Ptil * Ptil + 59.513639 * Ptil * Ptil * Ptil * Ptil
            //          + 9.6706522 * Ptil * Ptil * Ptil * Ptil * Ptil
            //          - 4.1136048 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //          - 1.0632301 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //          + 0.064400203 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //          + 0.035974396 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //          + 0.0025740066 * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil * Ptil
            //                * Ptil;
            // Tday  = Tnight;
            kt_hs = 0.0; //      Temperature forcing constant.
        }


        Tnight = 1.38877348e+03 + 2.79575848e+02 * Ptil - 2.13835822e+02 * pow(Ptil, 2)
                 + 2.10010475e+01 * pow(Ptil, 3) + 1.00938036e+02 * pow(Ptil, 4)
                 + 1.27972336e+01 * pow(Ptil, 5) - 1.39266925e+01 * pow(Ptil, 6)
                 - 3.70783272 * pow(Ptil, 7) + 5.22370269e-01 * pow(Ptil, 8)
                 + 3.20837882e-01 * pow(Ptil, 9) + 4.51831612e-02 * pow(Ptil, 10)
                 + 2.18195583e-03 * pow(Ptil, 11) + 3.98938097e-06 * pow(Ptil, 12);

        Tday = 2.15206036e+03 + 2.93485512e+01 * Ptil - 1.83318696e+02 * pow(Ptil, 2)
               + 4.63893130e+01 * pow(Ptil, 3) + 1.98116485e+01 * pow(Ptil, 4)
               - 2.85473177e+01 * pow(Ptil, 5) - 2.52726545 * pow(Ptil, 6)
               + 8.43627538 * pow(Ptil, 7) + 2.62945375 * pow(Ptil, 8)
               - 2.97098168e-01 * pow(Ptil, 9) - 2.86871487e-01 * pow(Ptil, 10)
               - 5.90629443e-02 * pow(Ptil, 11) - 5.38679474e-03 * pow(Ptil, 12)
               - 1.89972415e-04 * pow(Ptil, 13);

        if (pre < 100) {
            Tnight *= exp(0.1 * (log10(pre) - log10(100.0)));
            if (Tnight < 250)
                Tnight = 250;
            Tday *= exp(0.015 * (log10(pre) - log10(100.0)));
            if (Tday < 1000)
                Tday = 1000;
        }

        // Calculate Teq from eqn 26 (Heng, Menou, and Phillips 2011)
        if ((lon >= M_PI / 2) && (lon < 3 * M_PI / 2)) {
            Teq_hs = pow(Tnight * Tnight * Tnight * Tnight
                             + (Tday * Tday * Tday * Tday - Tnight * Tnight * Tnight * Tnight)
                                   * cos(lon - M_PI) * cos(lat),
                         0.25);
        }
        else {
            Teq_hs = Tnight;
        }

        //      Momentum dissipation constant.
        kv_hs = 0.0; //no boundary layer friction

        //      Update momenta
        for (int k = 0; k < 3; k++)
            Mh_d[id * 3 * nv + lev * 3 + k] =
                Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + kv_hs * time_step);

        //      Update temperature
        temperature_d[id * nv + lev] -= kt_hs * time_step * (temperature_d[id * nv + lev] - Teq_hs);
        // calculate heating rate as dP/dt
        // profx_dP_d[id * nv + lev] +=
        //     Rho_d[id * nv + lev] * Rd * (-kt_hs * (temperature_d[id * nv + lev] - Teq_hs));
    }
}
