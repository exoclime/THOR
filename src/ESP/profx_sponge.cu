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
// Description: Sponge Layer
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
#include <stdio.h>

#include "phy/profx_sponge.h"

__global__ void zonal_uv(double *M_d,
                         double *Rho_d,
                         int *   zonal_mean_tab_d,
                         double *lonlat_d,
                         int     num,
                         double *utmp,
                         double *vtmp,
                         int     max_count) {

    int    id  = blockIdx.x * blockDim.x + threadIdx.x;
    int    nv  = gridDim.y;
    int    lev = blockIdx.y;
    double u, v;

    if (id < num) {
        double lon;
        double rho = Rho_d[id * nv + lev];
        int    ind = zonal_mean_tab_d[id * 3];
        if (lonlat_d[id * 2] < 0)
            lon = lonlat_d[id * 2] + 2.0 * M_PI;
        else
            lon = lonlat_d[id * 2];
        u = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lon))
             + M_d[id * nv * 3 + lev * 3 + 1] * (cos(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * 0)
            / rho;

        v = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lonlat_d[id * 2 + 1]) * cos(lon))
             + M_d[id * nv * 3 + lev * 3 + 1] * (-sin(lonlat_d[id * 2 + 1]) * sin(lon))
             + M_d[id * nv * 3 + lev * 3 + 2] * (cos(lonlat_d[id * 2 + 1])))
            / rho;

        utmp[ind * nv * max_count + lev * max_count + zonal_mean_tab_d[id * 3 + 2]] =
            u / zonal_mean_tab_d[id * 3 + 1];
        vtmp[ind * nv * max_count + lev * max_count + zonal_mean_tab_d[id * 3 + 2]] =
            v / zonal_mean_tab_d[id * 3 + 1];
    }
}

__global__ void
zonal_w(double *W_d, double *Rho_d, int *zonal_mean_tab_d, int num, double *wtmp, int max_count) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double rho = Rho_d[id * nv + lev];
        int    ind = zonal_mean_tab_d[id * 3];

        wtmp[ind * nv * max_count + lev * max_count + zonal_mean_tab_d[id * 3 + 2]] =
            (W_d[id * nv + lev] / rho) / zonal_mean_tab_d[id * 3 + 1];
    }
}

__global__ void zonal_temp(double *pressure_d,
                           double *Rho_d,
                           double *Tbar_d,
                           int *   zonal_mean_tab_d,
                           double *lonlat_d,
                           int     num,
                           double *Ttmp,
                           double *Rd_d,
                           int     max_count) {

    //calculate zonal average temp for thermal sponge (aka Newtonian cooling)
    int    id  = blockIdx.x * blockDim.x + threadIdx.x;
    int    nv  = gridDim.y;
    int    lev = blockIdx.y;
    double temperature;

    if (id < num) {
        double rho = Rho_d[id * nv + lev];
        int    ind = zonal_mean_tab_d[id * 3];

        temperature = pressure_d[id * nv + lev] / Rd_d[id * nv + lev] / rho;

        Ttmp[ind * nv * max_count + lev * max_count + zonal_mean_tab_d[id * 3 + 2]] =
            temperature / zonal_mean_tab_d[id * 3 + 1];
    }
}

void print_vbar(double *vbar_h, int nlat, int nv) {
    int ind, lev = nv - 1;

    printf("\n");
    for (ind = 0; ind < nlat; ind++) {
        printf("%d, %g, %g, %g\n",
               ind,
               vbar_h[ind * nv * 3 + lev * 3 + 0],
               vbar_h[ind * nv * 3 + lev * 3 + 1],
               vbar_h[ind * nv * 3 + lev * 3 + 2]);
    }
}

__global__ void sponge_layer(double *M_d,
                             double *Rho_d,
                             double *W_d,
                             double *Wh_d,
                             double *pressure_d,
                             double *vbar_d,
                             double *Tbar_d,
                             int *   zonal_mean_tab_d,
                             double *lonlat_d,
                             double *Altitude_d,
                             double *Altitudeh_d,
                             double  Ruv,
                             double  Rw,
                             double  RT,
                             double  Rv_fac,
                             double  nsi,
                             bool    damp_uv_to_mean,
                             bool    damp_w_to_mean,
                             bool    implicit,
                             double  dt,
                             double *Rd_d,
                             int     nlat,
                             int     num,
                             int     nv,
                             bool    temp_sponge,
                             double *profx_dMh_d,
                             double *profx_dWh_d,
                             double *profx_dW_d,
                             double *profx_Qheat_d) {
    //
    //  Description:
    //
    //
    //
    //  Input: .
    //
    //  Output:
    //

    double n;
    double ns   = nsi;
    double ztop = Altitudeh_d[nv];
    double u, v, w;
    double vx, vy, vz;
    double rho;
    int    id = blockIdx.x * blockDim.x + threadIdx.x;

    if (id < num) {
        int    intlat  = zonal_mean_tab_d[id * 3];
        double des_lat = M_PI / nlat;
        double vbu, vbv, vbw;
        double temp_av;
        double lat1, lat2, lat3, lat;
        double intt, intl;
        lat1 = -M_PI / 2 + 0.5 * des_lat + des_lat * (intlat - 1);
        lat2 = -M_PI / 2 + 0.5 * des_lat + des_lat * (intlat);
        lat3 = -M_PI / 2 + 0.5 * des_lat + des_lat * (intlat + 1);
        lat  = lonlat_d[id * 2 + 1];

        for (int lev = 0; lev < nv; lev++) {
            rho = Rho_d[id * nv + lev];
            double lon;
            if (lonlat_d[id * 2] < 0)
                lon = lonlat_d[id * 2] + 2.0 * M_PI;
            else
                lon = lonlat_d[id * 2];

            if (damp_uv_to_mean) {
                // retrieve zonal averages and interpolate if necessary
                if (lat <= -M_PI / 2 + 0.5 * des_lat) {
                    //near poles, just use value in polar ring
                    vbu = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                    vbv = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                }
                else if (lat >= -M_PI / 2 + 0.5 * des_lat + des_lat * (nlat - 1)) {
                    //near poles, just use value in polar ring
                    vbu = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                    vbv = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                }
                else if (lat >= lat1 && lat < lat2) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat2) / (lat1 - lat2);
                    intl = (lat - lat1) / (lat2 - lat1);

                    vbu = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 0] * intt
                          + vbar_d[intlat * nv * 3 + lev * 3 + 0] * intl;
                    vbv = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 1] * intt
                          + vbar_d[intlat * nv * 3 + lev * 3 + 1] * intl;
                }
                else if (lat >= lat2 && lat < lat3) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat3) / (lat2 - lat3);
                    intl = (lat - lat2) / (lat3 - lat2);

                    vbu = vbar_d[(intlat)*nv * 3 + lev * 3 + 0] * intt
                          + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 0] * intl;
                    vbv = vbar_d[(intlat)*nv * 3 + lev * 3 + 1] * intt
                          + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 1] * intl;
                }
            }
            if (damp_w_to_mean) {
                // retrieve zonal averages and interpolate if necessary
                if (lat <= -M_PI / 2 + 0.5 * des_lat) {
                    //near poles, just use value in polar ring
                    vbw = vbar_d[intlat * nv * 3 + lev * 3 + 2];
                }
                else if (lat >= -M_PI / 2 + 0.5 * des_lat + des_lat * (nlat - 1)) {
                    //near poles, just use value in polar ring
                    vbw = vbar_d[intlat * nv * 3 + lev * 3 + 2];
                }
                else if (lat >= lat1 && lat < lat2) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat2) / (lat1 - lat2);
                    intl = (lat - lat1) / (lat2 - lat1);

                    vbw = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 2] * intt
                          + vbar_d[intlat * nv * 3 + lev * 3 + 2] * intl;
                }
                else if (lat >= lat2 && lat < lat3) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat3) / (lat2 - lat3);
                    intl = (lat - lat2) / (lat3 - lat2);

                    vbw = vbar_d[(intlat)*nv * 3 + lev * 3 + 2] * intt
                          + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 2] * intl;
                }
            }
            if (temp_sponge) {
                // retrieve zonal averages and interpolate if necessary
                if (lat <= -M_PI / 2 + 0.5 * des_lat) {
                    //near poles, just use value in polar ring
                    temp_av = Tbar_d[intlat * nv + lev];
                }
                else if (lat >= -M_PI / 2 + 0.5 * des_lat + des_lat * (nlat - 1)) {
                    //near poles, just use value in polar ring
                    temp_av = Tbar_d[intlat * nv + lev];
                }
                else if (lat >= lat1 && lat < lat2) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat2) / (lat1 - lat2);
                    intl = (lat - lat1) / (lat2 - lat1);

                    temp_av =
                        Tbar_d[(intlat - 1) * nv + lev] * intt + Tbar_d[intlat * nv + lev] * intl;
                }
                else if (lat >= lat2 && lat < lat3) {
                    // interpolate between nearest rings to smooth over latitude
                    intt = (lat - lat3) / (lat2 - lat3);
                    intl = (lat - lat2) / (lat3 - lat2);

                    temp_av =
                        Tbar_d[(intlat)*nv + lev] * intt + Tbar_d[(intlat + 1) * nv + lev] * intl;
                }
            }
            double kuv, kw;

            // calculate current winds
            u = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lon))
                 + M_d[id * nv * 3 + lev * 3 + 1] * (cos(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * 0)
                / rho;

            v = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lonlat_d[id * 2 + 1]) * cos(lon))
                 + M_d[id * nv * 3 + lev * 3 + 1] * (-sin(lonlat_d[id * 2 + 1]) * sin(lon))
                 + M_d[id * nv * 3 + lev * 3 + 2] * (cos(lonlat_d[id * 2 + 1])))
                / rho;

            w = W_d[id * nv + lev] / rho;

            // calculate strength of sponge
            n = Altitude_d[lev] / ztop;
            if (n >= ns) {
                kuv = Rv_fac * Ruv * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0); //
                kw  = Rv_fac * Rw * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0);  //
            }
            else {
                kuv = 0.0;
                kw  = 0.0;
            }

            //try implicit scheme
            if (implicit) {
                double unew, vnew, wnew;
                if (damp_uv_to_mean) {
                    unew = (u / dt + kuv * vbu) / (1.0 / dt + kuv);
                    vnew = (v / dt + kuv * vbv) / (1.0 / dt + kuv);
                }
                else {
                    unew = (u / dt) / (1.0 / dt + kuv);
                    vnew = (v / dt) / (1.0 / dt + kuv);
                }
                if (damp_w_to_mean) {
                    wnew = (w / dt + kw * vbw) / (1.0 / dt + kw);
                }
                else {
                    wnew = (w / dt) / (1.0 / dt + kw);
                }
                vx = unew * (-sin(lon)) + vnew * (-sin(lonlat_d[id * 2 + 1]) * cos(lon));
                vy = unew * (cos(lon)) + vnew * (-sin(lonlat_d[id * 2 + 1]) * sin(lon));
                vz = vnew * (cos(lonlat_d[id * 2 + 1]));
                M_d[id * nv * 3 + lev * 3 + 0] = vx * rho;
                M_d[id * nv * 3 + lev * 3 + 1] = vy * rho;
                M_d[id * nv * 3 + lev * 3 + 2] = vz * rho;
                W_d[id * nv + lev]             = wnew * rho;
            }
            else {
                double du, dv, dw;
                // calculate change in speed and convert to xyz
                if (damp_uv_to_mean) {
                    du = -kuv * (u - vbu);
                    dv = -kuv * (v - vbv);
                }
                else {
                    du = -kuv * u;
                    dv = -kuv * v;
                }

                vx = du * (-sin(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * cos(lon));

                vy = du * (cos(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * sin(lon));

                vz = dv * (cos(lonlat_d[id * 2 + 1]));

                profx_dMh_d[id * nv * 3 + lev * 3 + 0] +=
                    vx * rho; //drag rate... (not total change)
                profx_dMh_d[id * nv * 3 + lev * 3 + 1] += vy * rho;
                profx_dMh_d[id * nv * 3 + lev * 3 + 2] += vz * rho;

                // same for vertical speed
                if (damp_w_to_mean) {
                    dw = -kw * (w - vbw);
                }
                else {
                    dw = -kw * (w);
                }

                profx_dW_d[id * nv + lev] += dw * rho;
            }

            if (temp_sponge) {
                double temperature, dT, kvT;
                if (n >= ns) {
                    kvT =
                        Rv_fac * RT * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0); //
                }
                else {
                    kvT = 0.0;
                }
                temperature =
                    pressure_d[id * nv + lev] / Rd_d[id * nv + lev] / Rho_d[id * nv + lev];
                if (implicit) {
                    pressure_d[id * nv + lev] = Rho_d[id * nv + lev] * Rd_d[id * nv + lev]
                                                * (temperature / dt + kvT * temp_av)
                                                / (1.0 / dt + kvT);
                }
                else {
                    dT = -kvT * (temperature - temp_av);
                    profx_Qheat_d[id * nv + lev] += Rd_d[id * nv + lev] * Rho_d[id * nv + lev] * dT;
                }
            }
        }
        for (int lev = 1; lev < nv; lev++) {
            // interpolate vertical speed to boundaries
            double wl, wt;
            double xi, xim, xip;
            double intt, intl;

            xi  = Altitudeh_d[lev];
            xim = Altitude_d[lev - 1];
            xip = Altitude_d[lev];

            intt = (xi - xip) / (xim - xip);
            intl = (xi - xim) / (xip - xim);

            // update velocity directly
            if (implicit) {
                wl = W_d[id * nv + lev - 1];
                wt = W_d[id * nv + lev];

                Wh_d[id * (nv + 1) + lev] = wl * intt + wt * intl;
            }
            else {
                // send derivative to dynamical core
                wl = profx_dW_d[id * nv + lev - 1];
                wt = profx_dW_d[id * nv + lev];

                profx_dWh_d[id * (nv + 1) + lev] = wl * intt + wt * intl;
            }
        }
    }
}
