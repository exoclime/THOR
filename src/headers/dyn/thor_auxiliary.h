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
                                              double *Mh_d,
                                              double *h_d,
                                              double *hh_d,
                                              double *pt_d,
                                              double *pth_d,
                                              double *gtil_d,
                                              double *gtilh_d,
                                              double *Wh_d,
                                              double *W_d,
                                              double *epotential_d,
                                              double *epotentialh_d,
                                              double *ekinetic_d,
                                              double *ekinetich_d,
                                              double *Etotal_tau_d,
                                              double  P_Ref,
                                              double  Gravit,
                                              double *Cp_d,
                                              double *Rd_d,
                                              double *Altitude_d,
                                              double *Altitudeh_d,
                                              int     num,
                                              int     nv,
                                              bool    calcT,
                                              bool    energy_equation) {

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
    double temperature, rho, rhoh, rhol, rhop;
    double pressure, h, hl, pt, ptl, pl, pp, dpdz, ep, epl, ek, ekl;
    double gtilh, gtilht;
    double dz, dp_dz;
    double Cv;
    double CvoCp;

    double xi, xim, xip;
    double intt, intl, extr;
    double alt, altl, alth, altht;

    if (id < num) {
        for (int lev = 0; lev < nv + 1; lev++) {
            if (lev < nv) {
                Cv    = Cp_d[id * nv + lev] - Rd_d[id * nv + lev];
                CvoCp = Cv / Cp_d[id * nv + lev];
                if (lev > 0) {
                    pl   = pressure;
                    rhol = rho;
                    hl   = h;
                    if (energy_equation) {
                        epl = ep;
                        ekl = ek;
                    }
                    else {
                        ptl = pt;
                    }
                    altl = alt;
                }

                pressure = pressure_d[id * nv + lev];
                rho      = Rho_d[id * nv + lev];
                if (calcT) {
                    temperature                  = pressure / (rho * Rd_d[id * nv + lev]);
                    temperature_d[id * nv + lev] = temperature;
                }
                else {
                    temperature = temperature_d[id * nv + lev];
                }

                alt                = Altitude_d[lev];
                alth               = Altitudeh_d[lev];
                altht              = Altitudeh_d[lev + 1];
                h                  = Cv * temperature + pressure / rho; //why not just Cp*T??
                h_d[id * nv + lev] = h;

                if (energy_equation) {
                    //calculate energies needed later for total energy equation
                    ep                          = Gravit * Altitude_d[lev];
                    epotential_d[id * nv + lev] = ep;
                    ek                          = 0.5
                         * (Mh_d[id * nv * 3 + lev * 3 + 0] * Mh_d[id * nv * 3 + lev * 3 + 0]
                            + Mh_d[id * nv * 3 + lev * 3 + 1] * Mh_d[id * nv * 3 + lev * 3 + 1]
                            + Mh_d[id * nv * 3 + lev * 3 + 2] * Mh_d[id * nv * 3 + lev * 3 + 2]
                            + W_d[id * nv + lev] * W_d[id * nv + lev])
                         / (Rho_d[id * nv + lev] * Rho_d[id * nv + lev]);
                    ekinetic_d[id * nv + lev] = ek;
                    Etotal_tau_d[id * nv + lev] =
                        Rho_d[id * nv + lev] * (ep + ek) + Cv / Rd_d[id * nv + lev] * pressure;
                }
                else {
                    pt = (P_Ref / (Rd_d[id * nv + lev] * rho)) * pow(pressure / P_Ref, CvoCp);
                    pt_d[id * nv + lev] = pt;
                }
            }
            if (lev == 0) { // not sure about use of Rd in extrapolation
                extr = (-alt - Altitude_d[lev + 1]) / (Altitude_d[lev + 1] - alt);
                hh_d[id * (nv + 1) + 0] = h;

                pl   = pressure_d[id * nv + 1] - rho * Gravit * (-alt - Altitude_d[1]);
                rhol = Rho_d[id * nv + 1] + (Rho_d[id * nv + 1] - rho) * extr;
                rhoh = 0.5 * (rho + rhol);

                dz                         = 2.0 * Altitude_d[0];
                dpdz                       = (pressure_d[id * nv + 0] - pl) / dz;
                gtilh_d[id * (nv + 1) + 0] = -(1.0 / rhoh) * dpdz;
                gtilh                      = 0.0;
            }
            else if (lev == nv) {
                extr = (2 * Altitudeh_d[nv] - alt - Altitude_d[nv - 2])
                       / (Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                Cv                       = Cp_d[id * nv + lev - 1] - Rd_d[id * nv + lev - 1];
                CvoCp                    = Cv / Cp_d[id * nv + lev - 1];
                hh_d[id * (nv + 1) + nv] = h;

                pp = pressure_d[id * nv + nv - 2]
                     - rho * Gravit * (2 * Altitudeh_d[nv] - alt - Altitude_d[nv - 2]);
                if (pp < 0)
                    pp = 0;
                rhop = Rho_d[id * nv + nv - 2]
                       + (Rho_d[id * nv + nv - 1] - Rho_d[id * nv + nv - 2]) * extr;
                rhoh = 0.5 * (rho + rhop);

                if (energy_equation) {
                    epotentialh_d[id * (nv + 1) + nv] = Gravit * Altitudeh_d[nv];
                    epotentialh_d[id * (nv + 1) + 0]  = Gravit * Altitudeh_d[0];
                    ekinetich_d[id * (nv + 1) + nv] =
                        ekinetic_d[id * nv + nv - 2]
                        + (ekinetic_d[id * nv + nv - 1] - ekinetic_d[id * nv + nv - 2]) * extr;
                    ekinetich_d[id * (nv + 1) + 0] =
                        ekinetic_d[id * nv + 1]
                        + (ekinetic_d[id * nv + 1] - ekinetic_d[id * nv + 0])
                              * (Altitudeh_d[0] - Altitude_d[1]) / (Altitude_d[1] - Altitude_d[0]);
                }
                else {
                    pth_d[id * (nv + 1) + nv] =
                        pt_d[id * nv + nv - 2]
                        + (pt_d[id * nv + nv - 1] - pt_d[id * nv + nv - 2]) * extr;
                    pth_d[id * (nv + 1) + 0] = pt_d[id * nv + 1]
                                               + (pt_d[id * nv + 1] - pt_d[id * nv + 0])
                                                     * (Altitudeh_d[0] - Altitude_d[1])
                                                     / (Altitude_d[1] - Altitude_d[0]);
                }
                dz = 2.0 * (Altitudeh_d[nv] - Altitude_d[nv - 1]);
                pp = pressure_d[id * nv + nv - 2]
                     - Rho_d[id * nv + nv - 1] * Gravit
                           * (2 * Altitudeh_d[nv] - Altitude_d[nv - 1] - Altitude_d[nv - 2]);
                dpdz = (pp - pressure_d[id * nv + nv - 1]) / dz;

                gtilh_d[id * (nv + 1) + nv] = -(1.0 / rhoh) * dpdz;
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

                hh_d[id * (nv + 1) + lev] = hl * intt + h * intl;

                if (energy_equation) {
                    epotentialh_d[id * (nv + 1) + lev] = epl * intt + ep * intl;
                    ekinetich_d[id * (nv + 1) + lev]   = ekl * intt + ek * intl;
                }
                else {
                    pth_d[id * (nv + 1) + lev] = ptl * intt + pt * intl;
                }
                rhoh   = rhol * intt + rho * intl;
                gtilht = -(Wh_d[id * (nv + 1) + lev] / rhoh) * dp_dz;
                // since gtil_d (cell centers) goes only into slow modes,
                // we add velocity directly here

                xi  = alt;
                xim = alth;
                xip = altht;

                intt = (xi - xip) / (xim - xip);
                intl = (xi - xim) / (xip - xim);

                gtil_d[id * nv + lev - 1] = gtilh * intt + gtilht * intl;
                gtilh                     = gtilht;

                gtilh_d[id * (nv + 1) + lev] = -(1.0 / rhoh) * dp_dz;
                // since gtilh_d (interfaces) is in the coefficients of the vertical
                // equation, we don't want to include velocity
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

void add_fluxes_h(double *flux, double *areasT, int point_num, int nv) {
    int    i, lev;
    double flux_tot = 0.0, flux_max;

    for (lev = 0; lev < nv; lev++) {
        flux_tot = 0.0;
        flux_max = fabs(flux[0 * nv + lev]) * areasT[0];
        for (i = 0; i < point_num; i++) {
            flux_tot += flux[i * nv + lev] * areasT[i];
            if (fabs(flux[i * nv + lev] * areasT[i]) > flux_max) {
                flux_max = fabs(flux[i * nv + lev] * areasT[i]);
            }
        }
        printf("lev = %d, flux_tot = %e, flux_max = %e\n", lev, flux_tot, flux_max);
    }
}

__global__ void CalcMassThor(double *Mass_d,
                             double *GlobalMass_d,
                             double *Rho_d,
                             double  A,
                             double *Altitudeh_d,
                             double *lonlat_d,
                             double *areasT,
                             int     num,
                             bool    DeepModel) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        //calculate control volume
        double zup, zlow, Vol;
        zup  = Altitudeh_d[lev + 1] + A;
        zlow = Altitudeh_d[lev] + A;
        if (DeepModel) {
            Vol = areasT[id] / pow(A, 2) * (pow(zup, 3) - pow(zlow, 3)) / 3;
        }
        else {
            Vol = areasT[id] * (zup - zlow);
        }

        //mass in control volume = density*volume
        Mass_d[id * nv + lev] = (Rho_d[id * nv + lev]) * Vol;
    }
}
