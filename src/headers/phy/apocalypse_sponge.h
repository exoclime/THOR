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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////


__device__ double myatomicAdd(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int  old            = *address_as_ull, assumed;

    do {
        assumed = old;
        old     = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
        //  Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

__global__ void zonal_v(double *M_d,
                        double *W_d,
                        double *Rho_d,
                        double *vbar_d,
                        int *   zonal_mean_tab_d,
                        double *lonlat_d,
                        int     num) {

    int    id  = blockIdx.x * blockDim.x + threadIdx.x;
    int    nv  = gridDim.y;
    int    lev = blockIdx.y;
    double u, v;

    if (id < num) {
        double lon;
        double rho = Rho_d[id * nv + lev];
        int    ind = zonal_mean_tab_d[id * 2];
        if (lonlat_d[id * 2] < 0)
            lon = lonlat_d[id * 2] + 2.0 * M_PI;
        else
            lon = lonlat_d[id * 2];
        u = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lon)) + M_d[id * nv * 3 + lev * 3 + 1] * (cos(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * 0) / rho;

        v = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lonlat_d[id * 2 + 1]) * cos(lon)) + M_d[id * nv * 3 + lev * 3 + 1] * (-sin(lonlat_d[id * 2 + 1]) * sin(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * (cos(lonlat_d[id * 2 + 1]))) / rho;

        myatomicAdd(&(vbar_d[ind * nv * 3 + lev * 3 + 0]), u / zonal_mean_tab_d[id * 2 + 1]);
        myatomicAdd(&(vbar_d[ind * nv * 3 + lev * 3 + 1]), v / zonal_mean_tab_d[id * 2 + 1]);
        myatomicAdd(&(vbar_d[ind * nv * 3 + lev * 3 + 2]), (W_d[id * nv + lev] / rho) / zonal_mean_tab_d[id * 2 + 1]);
    }
}

__global__ void sponge_layer(double *M_d,
                             double *Rho_d,
                             double *W_d,
                             double *Wh_d,
                             double *vbar_d,
                             int *   zonal_mean_tab_d,
                             double *lonlat_d,
                             double *Altitude_d,
                             double *Altitudeh_d,
                             double  Rv,
                             double  nsi,
                             double  dt,
                             int     nlat,
                             int     num,
                             int     nv) {
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
        int    intlat  = zonal_mean_tab_d[id * 2];
        double des_lat = M_PI / nlat;
        double vbu, vbv, vbw;
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

            if (lat <= -M_PI / 2 + 0.5 * des_lat) {
                vbu = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                vbv = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                vbw = vbar_d[intlat * nv * 3 + lev * 3 + 2];
            }
            else if (lat >= -M_PI / 2 + 0.5 * des_lat + des_lat * (nlat - 1)) {
                vbu = vbar_d[intlat * nv * 3 + lev * 3 + 0];
                vbv = vbar_d[intlat * nv * 3 + lev * 3 + 1];
                vbw = vbar_d[intlat * nv * 3 + lev * 3 + 2];
            }
            else if (lat >= lat1 && lat < lat2) {

                intt = (lat - lat2) / (lat1 - lat2);
                intl = (lat - lat1) / (lat2 - lat1);

                vbu = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 0] * intt + vbar_d[intlat * nv * 3 + lev * 3 + 0] * intl;
                vbv = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 1] * intt + vbar_d[intlat * nv * 3 + lev * 3 + 1] * intl;
                vbw = vbar_d[(intlat - 1) * nv * 3 + lev * 3 + 2] * intt + vbar_d[intlat * nv * 3 + lev * 3 + 2] * intl;
            }
            else if (lat >= lat2 && lat < lat3) {

                intt = (lat - lat3) / (lat2 - lat3);
                intl = (lat - lat2) / (lat3 - lat2);

                vbu = vbar_d[(intlat)*nv * 3 + lev * 3 + 0] * intt + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 0] * intl;
                vbv = vbar_d[(intlat)*nv * 3 + lev * 3 + 1] * intt + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 1] * intl;
                vbw = vbar_d[(intlat)*nv * 3 + lev * 3 + 2] * intt + vbar_d[(intlat + 1) * nv * 3 + lev * 3 + 2] * intl;
            }
            double kv;

            u = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lon)) + M_d[id * nv * 3 + lev * 3 + 1] * (cos(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * 0) / rho;

            v = (M_d[id * nv * 3 + lev * 3 + 0] * (-sin(lonlat_d[id * 2 + 1]) * cos(lon)) + M_d[id * nv * 3 + lev * 3 + 1] * (-sin(lonlat_d[id * 2 + 1]) * sin(lon)) + M_d[id * nv * 3 + lev * 3 + 2] * (cos(lonlat_d[id * 2 + 1]))) / rho;
            ;

            w = W_d[id * nv + lev] / rho;

            n = Altitude_d[lev] / ztop;
            if (n >= ns) {
                kv = Rv * pow(sin(0.5 * M_PI * (n - ns) * (1.0) / (1.0 - ns)), 2.0);
            }
            else {
                kv = 0.0;
            }

            du = -kv * (u - vbu) * dt;
            dv = -kv * (v - vbv) * dt;

            vx = du * (-sin(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * cos(lon));

            vy = du * (cos(lon)) + dv * (-sin(lonlat_d[id * 2 + 1]) * sin(lon));

            vz = dv * (cos(lonlat_d[id * 2 + 1]));

            M_d[id * nv * 3 + lev * 3 + 0] += vx * rho;
            M_d[id * nv * 3 + lev * 3 + 1] += vy * rho;
            M_d[id * nv * 3 + lev * 3 + 2] += vz * rho;

            dw = -kv * (w - vbw) * dt;
            W_d[id * nv + lev] += dw * rho;
        }
        for (int lev = 1; lev < nv; lev++) {

            double wl, wt;
            double xi, xim, xip;
            double intt, intl;

            xi  = Altitudeh_d[lev];
            xim = Altitude_d[lev - 1];
            xip = Altitude_d[lev];

            intt = (xi - xip) / (xim - xip);
            intl = (xi - xim) / (xip - xim);

            wl = W_d[id * nv + lev - 1];
            wt = W_d[id * nv + lev];

            Wh_d[id * (nv + 1) + lev] = wl * intt + wt * intl;
        }
    }
}
