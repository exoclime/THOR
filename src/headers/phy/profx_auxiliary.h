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
// Description: Computes temperature and pressure.
//
//
// Method: We use the assumption that the atmosphere can be treated as an ideal gas.
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

__global__ void Compute_temperature_only(double *temperature_d,
                                         double *pressure_d,
                                         double *Rho_d,
                                         double *Rd_d,
                                         int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    // Computes absolute temperature
    if (id < num)
        temperature_d[id * nv + lev] =
            pressure_d[id * nv + lev] / (Rd_d[id * nv + lev] * Rho_d[id * nv + lev]);
}

__global__ void Compute_temperature(double *temperature_d,
                                    double *pt_d,
                                    double *pressure_d,
                                    double *Rho_d,
                                    double  P_Ref,
                                    double *Rd_d,
                                    double *Cp_d,
                                    int     num,
                                    bool    calcT) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    double Cv;
    double CvoCp;

    // Computes absolute and potential temperature
    if (id < num) {
        Cv    = Cp_d[id * nv + lev] - Rd_d[id * nv + lev];
        CvoCp = Cv / Cp_d[id * nv + lev];
        if (calcT) {
            temperature_d[id * nv + lev] =
                pressure_d[id * nv + lev] / (Rd_d[id * nv + lev] * Rho_d[id * nv + lev]);
        }
        pt_d[id * nv + lev] = (P_Ref / (Rd_d[id * nv + lev] * Rho_d[id * nv + lev]))
                              * pow(pressure_d[id * nv + lev] / P_Ref, CvoCp);
    }
}

__global__ void
Compute_pressure(double *pressure_d, double *temperature_d, double *Rho_d, double *Rd_d, int num) {


    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    // Computes absolute pressure
    if (id < num) {
        pressure_d[id * nv + lev] =
            temperature_d[id * nv + lev] * Rd_d[id * nv + lev] * Rho_d[id * nv + lev];
    }
}

__global__ void
Recompute_W(double *W_d, double *Wh_d, double *Altitude_d, double *Altitudeh_d, int num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double xi, xim, xip, whl, wht, intt, intl;
        whl = Wh_d[id * (nv + 1) + lev];
        wht = Wh_d[id * (nv + 1) + lev + 1];

        xi  = Altitude_d[lev];
        xim = Altitudeh_d[lev];
        xip = Altitudeh_d[lev + 1];

        intt = (xi - xip) / (xim - xip);
        intl = (xi - xim) / (xip - xim);

        W_d[id * nv + lev] = whl * intt + wht * intl;
    }
}

__global__ void isnan_check(double *array, int width, int height, bool *check) {

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

// __global__ void isnan_loop(double *t, int num, int nv) {
//   int i;
//
//   for (i=0;i<num;i++) {
//     if (isnan(t[i*nv])) {
//       printf("id = %d, lev = %d, t = %f\n",i,0,t[i*nv]);
//     }
//   }
//

__global__ void apply_heating(double *temperature_d,
                              double *profx_Qheat_d,
                              double *Rho_d,
                              double *Cp_d,
                              double *Rd_d,
                              double  timestep,
                              int     num) {
    //updates temperature from Qheat in the case that gcm_off = true
    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        temperature_d[id * nv + lev] += 1.0 / (Cp_d[id * nv + lev] - Rd_d[id * nv + lev])
                                        * profx_Qheat_d[id * nv + lev] / Rho_d[id * nv + lev]
                                        * timestep;
        if (temperature_d[id * nv + lev] < 0)
            temperature_d[id * nv + lev] = 0.0;
        if (isnan(temperature_d[id * nv + lev])) {
            printf("check\n");
        }
    }
}
