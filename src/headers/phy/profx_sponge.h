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


__global__ void zonal_v(double *M_d,
                        double *W_d,
                        double *Rho_d,
                        double *vbar_d,
                        int *   zonal_mean_tab_d,
                        double *lonlat_d,
                        int     num,
                        double *utmp,
                        double *vtmp,
                        double *wtmp,
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
                           double  Rd,
                           int     max_count) {

    //calculate zonal average temp for thermal sponge (aka Newtonian cooling)
    int    id  = blockIdx.x * blockDim.x + threadIdx.x;
    int    nv  = gridDim.y;
    int    lev = blockIdx.y;
    double temperature;

    if (id < num) {
        double rho = Rho_d[id * nv + lev];
        int    ind = zonal_mean_tab_d[id * 3];

        temperature = pressure_d[id * nv + lev] / Rd / rho;

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
                             double  Rv,
                             double  RvT,
                             double  Rv_fac,
                             double  nsi,
                             double  dt,
                             double  Rd,
                             int     nlat,
                             int     num,
                             int     nv,
                             bool    temp_sponge,
                             double *profx_dMh_d,
                             double *profx_dWh_d,
                             double *profx_dW_d) {
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
    double du, dv, dw;
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

            // retrieve zonal averages and interpolate if necessary
            if (lat <= -M_PI / 2 + 0.5 * des_lat) {
                //near poles, just use value in polar ring
                vbu     = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                vbv     = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                vbw     = vbar_d[intlat * nv * 3 + lev * 3 + 2];
                temp_av = Tbar_d[intlat * nv + lev];
            }
            else if (lat >= -M_PI / 2 + 0.5 * des_lat + des_lat * (nlat - 1)) {
                //near poles, just use value in polar ring
                vbu     = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                vbv     = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                vbw     = vbar_d[intlat * nv * 3 + lev * 3 + 2];
                temp_av = Tbar_d[intlat * nv + lev];
            }
            else if (lat >= lat1 && lat < lat2) {
                // interpolate between nearest rings to smooth over latitude
                intt = (lat - lat2) / (lat1 - lat2);
                intl = (lat - lat1) / (lat2 - lat1);

                vbu = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 0] * intt
                      + vbar_d[intlat * nv * 3 + lev * 3 + 0] * intl;
                vbv = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 1] * intt
                      + vbar_d[intlat * nv * 3 + lev * 3 + 1] * intl;
                vbw = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 2] * intt
                      + vbar_d[intlat * nv * 3 + lev * 3 + 2] * intl;
                temp_av = Tbar_d[(intlat - 1) * nv + lev] * intt + Tbar_d[intlat * nv + lev] * intl;
            }
            else if (lat >= lat2 && lat < lat3) {
                // interpolate between nearest rings to smooth over latitude
                intt = (lat - lat3) / (lat2 - lat3);
                intl = (lat - lat2) / (lat3 - lat2);

                vbu = vbar_d[(intlat)*nv * 3 + lev * 3 + 0] * intt
                      + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 0] * intl;
                vbv = vbar_d[(intlat)*nv * 3 + lev * 3 + 1] * intt
                      + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 1] * intl;
                vbw = vbar_d[(intlat)*nv * 3 + lev * 3 + 2] * intt
                      + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 2] * intl;
                temp_av = Tbar_d[(intlat)*nv + lev] * intt + Tbar_d[(intlat + 1) * nv + lev] * intl;
            }
            double kv;

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
                kv = Rv_fac * Rv * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0); //
            }
            else {
                kv = 0.0;
            }

            // calculate change in speed and convert to xyz
            du = -kv * (u - vbu); // * dt;
            dv = -kv * (v - vbv); // * dt;

            vx = du * (-sin(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * cos(lon));

            vy = du * (cos(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * sin(lon));

            vz = dv * (cos(lonlat_d[id * 2 + 1]));

            // update the momentum variables
            // M_d[id * nv * 3 + lev * 3 + 0] += vx * rho * dt;
            // M_d[id * nv * 3 + lev * 3 + 1] += vy * rho * dt;
            // M_d[id * nv * 3 + lev * 3 + 2] += vz * rho * dt;
            profx_dMh_d[id * nv * 3 + lev * 3 + 0] += vx * rho; //drag rate... (not total change)
            profx_dMh_d[id * nv * 3 + lev * 3 + 1] += vy * rho;
            profx_dMh_d[id * nv * 3 + lev * 3 + 2] += vz * rho;

            // same for vertical speed
            dw = -kv * (w - vbw); // * dt;
            // dw = -kv * (w); // * dt;

            profx_dW_d[id * nv + lev] += dw * rho;
            // dw = -kv * (w)*dt; //what happens if i damp to zero ??
            // W_d[id * nv + lev] += dw * rho * dt;

            if (temp_sponge) {
                double temperature, dT, kvT;
                if (n >= ns) {
                    kvT =
                        Rv_fac * RvT * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0); //
                }
                else {
                    kvT = 0.0;
                }
                temperature = pressure_d[id * nv + lev] / Rd / Rho_d[id * nv + lev];
                dT          = -kvT * (temperature - temp_av) * dt;
                pressure_d[id * nv + lev] += Rd * Rho_d[id * nv + lev] * dT;
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

            // wl = W_d[id * nv + lev - 1];
            // wt = W_d[id * nv + lev];
            //
            // Wh_d[id * (nv + 1) + lev] = wl * intt + wt * intl;

            wl = profx_dW_d[id * nv + lev - 1];
            wt = profx_dW_d[id * nv + lev];

            profx_dWh_d[id * (nv + 1) + lev] = wl * intt + wt * intl;
        }
    }
}
