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
// Description: Contains the subroutines that compute: - Temperature, enthalpy, potential temperature and effective gravity.
//
//
// Method: -
//
// Known limitations: None
//
// Known issues: None
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

#include "dyn/phy_modules_device.h"

extern __constant__     device_RK_array
                        dynamical_core_phy_modules_arrays[NUM_PHY_MODULES_DYN_CORE_ARRAYS];
extern __constant__ int num_dynamical_arrays[1];


__global__ void Compute_Temperature_H_Pt_Geff(double *temperature_d,
                                              double *pressure_d,
                                              double *Rho_d,
                                              double *h_d,
                                              double *hh_d,
                                              double *pt_d,
                                              double *pth_d,
                                              double *gtil_d,
                                              double *gtilh_d,
                                              double *Wh_d,
                                              double  P_Ref,
                                              double  Gravit,
                                              double  Cp,
                                              double  Rd,
                                              double *Altitude_d,
                                              double *Altitudeh_d,
                                              int     num,
                                              int     nv) {

    //
    //  Description: Computes temperature, internal energy, potential temperature and effective gravity.
    //
    //
    //  Input: density, vertical momentum, heat capacity, gravity, reference pressure, gas constant and altitudes.
    //
    //  Output: temperature, enthalpy, potential temperature and effective gravity.
    //
    //  Thread index.
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    //  Local variables.
    double temperature, rho, rhoh, rhol;
    double pressure, h, hl, pt, ptl, pl, pp, dpdz;
    double gtilh, gtilht;
    double dz, dp_dz;
    double Cv    = Cp - Rd;
    double CvoCp = Cv / Cp;

    double xi, xim, xip;
    double intt, intl;
    double alt, altl, alth, altht;

    if (id < num) {
        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev < nv) {
                if (lev > 0) {
                    pl   = pressure;
                    rhol = rho;
                    hl   = h;
                    ptl  = pt;
                    altl = alt;
                }

                pressure    = pressure_d[id * nv + lev];
                rho         = Rho_d[id * nv + lev];
                temperature = pressure / (Rd * rho);

                alt   = Altitude_d[lev];
                alth  = Altitudeh_d[lev];
                altht = Altitudeh_d[lev + 1];
                h     = Cv * temperature + pressure / rho;
                pt    = (P_Ref / (Rd * rho)) * pow(pressure / P_Ref, CvoCp);

                temperature_d[id * nv + lev] = temperature;
                h_d[id * nv + lev]           = h;
                pt_d[id * nv + lev]          = pt;
            }
            if (lev == 0) {
                hh_d[id * (nv + 1) + 0] = h;

                pl   = pressure_d[id * nv + 1] - rho * Gravit * (-alt - Altitude_d[1]);
                rhoh = 0.5 * (pressure + pl) / (Rd * temperature);
                pth_d[id * (nv + 1)] =
                    (P_Ref / (Rd * rhoh)) * pow(0.5 * (pressure + pl) / P_Ref, CvoCp);

                gtilh_d[id * (nv + 1)] = 0.0;
                gtilh                  = 0.0;
            }
            else if (lev == nv) {
                hh_d[id * (nv + 1) + nv] = h;

                pp = pressure_d[id * nv + nv - 2]
                     - rho * Gravit * (2 * Altitudeh_d[nv] - alt - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0;
                rhoh = 0.5 * (pressure + pp) / (Rd * temperature);
                pth_d[id * (nv + 1) + nv] =
                    (P_Ref / (Rd * rhoh)) * pow(0.5 * (pressure + pp) / P_Ref, CvoCp);

                gtilh_d[id * (nv + 1) + nv] = 0.0;
                gtilht                      = 0.0;

                xi  = alt;
                xim = alth;
                xip = altht;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                gtil_d[id * nv + nv - 1] = gtilh * intt + gtilht * intl;
            }
            else {

                dz    = alt - altl;
                dp_dz = (pressure - pl) / dz;

                xi  = alth;
                xim = altl;
                xip = alt;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                hh_d[id * (nv + 1) + lev]  = hl * intt + h * intl;
                pth_d[id * (nv + 1) + lev] = ptl * intt + pt * intl;

                rhoh   = rhol * intt + rho * intl;
                gtilht = -(Wh_d[id * (nv + 1) + lev] / rhoh) * dp_dz;

                gtilh_d[id * (nv + 1) + lev] = gtilht;

                xi  = alt;
                xim = alth;
                xip = altht;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);


                gtil_d[id * nv + lev - 1] = gtilh * intt + gtilht * intl;
                gtilh                     = gtilht;
            }
        }


        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev == 0) {
                dz = 2.0 * Altitude_d[0];
                pl = pressure_d[id * nv + 1]
                     - Rho_d[id * nv + 0] * Gravit * (-Altitude_d[0] - Altitude_d[1]);
                dpdz = (pressure_d[id * nv + 0] - pl) / dz;
                rhoh = 0.5 * (pressure_d[id * nv + 0] + pl) / (Rd * temperature_d[id * nv + 0]);

                gtilh_d[id * (nv + 1) + 0] = -(1.0 / rhoh) * dpdz;
            }
            else if (lev == nv) {

                dz = 2.0 * (Altitudeh_d[nv] - Altitude_d[nv - 1]);
                pp = pressure_d[id * nv + nv - 2]
                     - Rho_d[id * nv + nv - 1] * Gravit
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                dpdz = (pp - pressure_d[id * nv + nv - 1]) / dz;

                rhoh = 0.5 * (pressure_d[id * nv + nv - 1] + pp)
                       / (Rd * temperature_d[id * nv + nv - 1]);

                gtilh_d[id * (nv + 1) + nv] = -(1.0 / rhoh) * dpdz;
            }
            else {
                dz   = Altitude_d[lev] - Altitude_d[lev - 1];
                dpdz = (pressure_d[id * nv + lev] - pressure_d[id * nv + lev - 1]) / dz;

                xi  = Altitudeh_d[lev];
                xim = Altitude_d[lev - 1];
                xip = Altitude_d[lev];

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                rhoh = Rho_d[id * nv + lev - 1] * intt + Rho_d[id * nv + lev] * intl;

                gtilh_d[id * (nv + 1) + lev] = -(1.0 / rhoh) * dpdz;
            }
        }
    }
}

__global__ void UpdateRK(double *M_d,
                         double *Mk_d,
                         double *Mi_d,
                         double *Wh_d,
                         double *Whk_d,
                         double *Whi_d,
                         double *W_d,
                         double *Rho_d,
                         double *Rhok_d,
                         double *Rhoi_d,
                         double *pressure_d,
                         double *pressurek_d,
                         double *pressurei_d,
                         double *func_r_d,
                         double *Altitude_d,
                         double *Altitudeh_d,
                         int     num,
                         int     nv) {

    //
    //  Description: Updates deviations.
    //
    //  Input: initial arrays and arrays from the time integration step.
    //
    //  Output: Deviations for each diagnostic.
    //

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    double aux;
    double Mx, My, Mz;
    double fx, fy, fz;

    double xi, xim, xip;
    double althl, altht;

    double intt, intl;
    double wht, whl;

    if (id < num) {

        fx = func_r_d[id * 3 + 0];
        fy = func_r_d[id * 3 + 1];
        fz = func_r_d[id * 3 + 2];

        for (int lev = 0; lev < nv; lev++) {
            // Mh
            Mx = Mi_d[id * 3 * nv + lev * 3 + 0] - Mk_d[id * 3 * nv + lev * 3 + 0];
            My = Mi_d[id * 3 * nv + lev * 3 + 1] - Mk_d[id * 3 * nv + lev * 3 + 1];
            Mz = Mi_d[id * 3 * nv + lev * 3 + 2] - Mk_d[id * 3 * nv + lev * 3 + 2];

            // Radial component.
            aux = Mx * fx + My * fy + Mz * fz;

            // Subtract the radial component of the vector, if any.
            M_d[id * 3 * nv + lev * 3 + 0] = Mx - aux * fx;
            M_d[id * 3 * nv + lev * 3 + 1] = My - aux * fy;
            M_d[id * 3 * nv + lev * 3 + 2] = Mz - aux * fz;

            // Rho
            Rho_d[id * nv + lev] = Rhoi_d[id * nv + lev] - Rhok_d[id * nv + lev];

            // Pressure
            pressure_d[id * nv + lev] = pressurei_d[id * nv + lev] - pressurek_d[id * nv + lev];

            // phy modules arrays from array defs
            for (int i = 0; i < num_dynamical_arrays[0]; i++) {
                int     num_dim = dynamical_core_phy_modules_arrays[i].dimensions;
                double *array   = dynamical_core_phy_modules_arrays[i].array_d;
                double *arrayi  = dynamical_core_phy_modules_arrays[i].arrayi_d;
                double *arrayk  = dynamical_core_phy_modules_arrays[i].arrayk_d;
                for (int itr = 0; itr < num_dim; itr++) {
                    array[(id * nv + lev) * num_dim + itr] =
                        arrayi[(id * nv + lev) * num_dim + itr]
                        - arrayk[(id * nv + lev) * num_dim + itr];
                }
            }
        }
        // Wh
        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev == 0) {
                Wh_d[id * (nv + 1) + 0] = 0.0;
                altht                   = Altitudeh_d[lev];
                wht                     = 0.0;
            }
            else {
                whl = wht;
                if (lev == nv)
                    wht = 0.0;
                else
                    wht = Whi_d[id * (nv + 1) + lev] - Whk_d[id * (nv + 1) + lev];
                Wh_d[id * (nv + 1) + lev] = wht;

                althl = altht;
                altht = Altitudeh_d[lev];

                xi  = Altitude_d[lev - 1];
                xim = althl;
                xip = altht;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                W_d[id * nv + lev - 1] = whl * intt + wht * intl;
            }
        }
    }
}

__global__ void UpdateRK2(double *M_d,
                          double *Mk_d,
                          double *Wh_d,
                          double *Whk_d,
                          double *Wk_d,
                          double *Rho_d,
                          double *Rhok_d,
                          double *pressure_d,
                          double *pressurek_d,
                          double *func_r_d,
                          double *Altitude_d,
                          double *Altitudeh_d,
                          int     num,
                          int     nv) {

    //
    //  Description: Update arrays from the large step.
    //
    //  Input: deviations.
    //
    //  Output: diagnostics from the large time step.
    //

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    double aux;
    double Mx, My, Mz;
    double fx, fy, fz;

    double xi, xim, xip;
    double althl, altht;

    double intt, intl;
    double wht, whl;


    if (id < num) {

        fx = func_r_d[id * 3 + 0];
        fy = func_r_d[id * 3 + 1];
        fz = func_r_d[id * 3 + 2];

        for (int lev = 0; lev < nv; lev++) {
            // Mh
            Mx = M_d[id * 3 * nv + lev * 3 + 0] + Mk_d[id * 3 * nv + lev * 3 + 0];
            My = M_d[id * 3 * nv + lev * 3 + 1] + Mk_d[id * 3 * nv + lev * 3 + 1];
            Mz = M_d[id * 3 * nv + lev * 3 + 2] + Mk_d[id * 3 * nv + lev * 3 + 2];

            // Radial component.
            aux = Mx * fx + My * fy + Mz * fz;

            // Subtract the radial component of the vector, if any.
            Mk_d[id * 3 * nv + lev * 3 + 0] = Mx - aux * fx;
            Mk_d[id * 3 * nv + lev * 3 + 1] = My - aux * fy;
            Mk_d[id * 3 * nv + lev * 3 + 2] = Mz - aux * fz;

            // Rho
            Rhok_d[id * nv + lev] += Rho_d[id * nv + lev];

            // Pressure
            pressurek_d[id * nv + lev] += pressure_d[id * nv + lev];


            int sz = num_dynamical_arrays[0];
            // phy modules arrays from array defs
            for (int i = 0; i < sz; i++) {
                int     num_dim = dynamical_core_phy_modules_arrays[i].dimensions;
                double *array   = dynamical_core_phy_modules_arrays[i].array_d;
                double *arrayk  = dynamical_core_phy_modules_arrays[i].arrayk_d;

                for (int itr = 0; itr < num_dim; itr++) {
                    arrayk[(id * nv + lev) * num_dim + itr] +=
                        array[(id * nv + lev) * num_dim + itr];
                }
            }
        }
        // Wh
        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev == 0 || lev == nv)
                Whk_d[id * (nv + 1) + lev] = 0.0;
            else
                Whk_d[id * (nv + 1) + lev] += Wh_d[id * (nv + 1) + lev];
        }
        for (int lev = 0; lev < nv; lev++) {
            if (lev == 0) {
                althl = 0;
                altht = Altitudeh_d[1];
                whl   = 0.0;
                wht   = Whk_d[id * (nv + 1) + 1];
            }
            xi  = Altitude_d[lev];
            xim = althl;
            xip = altht;

            intt = (xi - xip) / (xim - xip);
            intl = (xi - xim) / (xip - xim);

            Wk_d[id * nv + lev] = whl * intt + wht * intl;
            if (lev < nv - 1) {
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                whl   = wht;
                wht   = Whk_d[id * (nv + 1) + lev + 2];
            }
        }
    }
}

__global__ void isnan_check_thor(double *array, int width, int height, bool *check) {

    //
    //  Description: Check for nan in the temperature array.
    //

    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    while (idx < width) {
        for (int i = 0; i < height; i++)
            if (isnan(array[(i * width) + idx])) {
                *check = true;
            }
        idx += gridDim.x + blockDim.x;
    }
}
