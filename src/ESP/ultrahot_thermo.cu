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
// ultrahot jupiter thermodynamics
//
//
// Known limitations: - Runs in a single GPU.
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

#include "esp.h"
#include "phy/ultrahot_thermo.h"


bool ESP::read_in_gibbs_H(int GibbsN) {
    FILE *infile;
    int   nlines = GibbsN; //lines in GibbsH file

    infile  = fopen(Gibbs_Filename, "r");
    GibbsT  = (double *)malloc(nlines * sizeof(double));
    GibbsdG = (double *)malloc(nlines * sizeof(double));

    cudaMalloc((void **)&GibbsT_d, nlines * sizeof(double));
    cudaMalloc((void **)&GibbsdG_d, nlines * sizeof(double));

    if (!path_exists(Gibbs_Filename)) {
        log::printf("\nGibbs H input file %s does not exist.\n", Gibbs_Filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < nlines; i++) {
        if (fscanf(infile, "%lf %lf", &GibbsT[i], &GibbsdG[i]) != 2) {
            log::printf("error parsing gibbs H file %s.\n", Gibbs_Filename);
            fclose(infile);
            return false;
        }
    }

    fclose(infile);
    return true;
}

__host__ __device__ double
         linear_interpolation(double x1, double x2, double y1, double y2, double x) {
    // calculates value of y at x on interval (x1, x2), given y1(x1), y2(x2)
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

__host__ __device__ int locate_min_i(double *array_d, int N, double val) {

    int id = -1;
    if (val >= array_d[N - 1]) {
        id = N - 2;
    }
    else if (val < array_d[0]) {
        id = 0;
    }
    else {
        for (int j = 1; j < N; j++) {
            if (val >= array_d[j - 1] && val < array_d[j]) {
                id = j - 1;
                break;
            }
        }
    }
    return id;
}

__host__ __device__ int locate_max_i(double *array_d, int N, double val) {

    int id = -1;
    if (val >= array_d[N - 1]) {
        id = N - 1;
    }
    else if (val < array_d[0]) {
        id = 1;
    }
    else {
        for (int j = 1; j < N; j++) {
            if (val >= array_d[j - 1] && val < array_d[j]) {
                id = j;
                break;
            }
        }
    }
    return id;
}

__host__ double guillot_T(double pressure, double mu, double Teq, double P_Ref, double Gravit) {
    // Guillot profile (temperature is a function of optical depth/pressure)
    // Based on Equation 27 of Guillot 2010 and Equation 41 of Heng+ 2011b

    //somehow need to set this up to use user-set RT params
    double Tint    = 800, T4;
    double tau_lw0 = 1998;
    double fl      = 0.5;
    double tau0    = fl * tau_lw0;
    double k0      = tau0 * Gravit / P_Ref;
    double eps     = 1.0 / fl;
    double tau_sw  = tau0 * 0.5;
    double ksw     = tau_sw * Gravit / P_Ref;
    double x0      = P_Ref / Gravit;
    double gamma0  = ksw / k0;
    //  double mu      = 0.5;

    double klw    = k0 * (1.0 + 2 * (eps - 1) * pressure / P_Ref);
    double x      = pressure / Gravit;
    double tau    = k0 * x;
    double tau_lw = tau0 * pressure / P_Ref + (eps - 1) * tau0 * pow(pressure / P_Ref, 2);
    double gamma  = ksw / klw;

    T4 = 0.75 * pow(Tint, 4) * (2. / 3 + tau_lw)
         + 3 * mu * pow(Teq, 4)
               * (2. / 3 + gamma / (3 * mu) * exp(-gamma0 * tau / mu)
                  + mu / gamma0 * (1.0 - exp(-gamma0 * tau / mu))
                  + 2 * (eps - 1) * mu / (gamma0 * ksw * x0)
                        * (mu - (mu + gamma0 * tau) * exp(-gamma0 * tau / mu)));
    return pow(T4, 0.25);
}

__host__ __device__ double chi_H_equilibrium(double *GibbsT,
                                             double *GibbsdG,
                                             int     GibbsN,
                                             double  temperature,
                                             double  pressure) {
    int imin, imax;
    imin = locate_min_i(GibbsT, GibbsN, temperature);
    imax = locate_max_i(GibbsT, GibbsN, temperature);
    double dGnew =
        linear_interpolation(GibbsT[imin], GibbsT[imax], GibbsdG[imin], GibbsdG[imax], temperature);
    double kprime = pressure / 1e5 * exp(2000 * dGnew / RUNIV / temperature);
    // printf("k = %f ", kprime);
    double n_H;
    if (isinf(kprime)) { //low temperature, no atomic H
        n_H = 0.0;
    }
    else if (kprime == 0) { //high temperature, all atomic H
        n_H = 1.0;
    }
    else {
        n_H = (-1.0 + sqrt(1 + 8 * kprime)) / (4 * kprime);
    }
    double X_H   = 2 * n_H / (n_H + 1);
    double mmean = 2 * massH - massH * X_H;
    return massH / mmean
           * X_H; //hmm, in this simple case, chi_H = n_H... wonder if there's any reason to generalize
}

__device__ __host__ double Rd_from_chi_H(double chi_H) {
    return kBoltz * chi_H / massH + kBoltz * (1 - chi_H) / (2 * massH);
}

__host__ double heat_capacity_H2(double temperature) {
    // heat capacity at constant pressure of molecular hydrogen from chase 1998 (nist database)
    double A, B, C, D, E, t = temperature / 1000.0;

    if (temperature < 1000) {
        A = 33.066178;
        B = -11.363417;
        C = 11.432816;
        D = -2.772874;
        E = -0.158558;
    }
    else if (temperature >= 1000 && temperature < 2500) {
        A = 18.563083;
        B = 12.257357;
        C = -2.859786;
        D = 0.268238;
        E = 1.977990;
    }
    else {
        A = 43.423560;
        B = -4.293079;
        C = 1.272428;
        D = -0.096876;
        E = -20.533862;
    }
    return (A + B * t + C * pow(t, 2) + D * pow(t, 3) + E * pow(t, -2)) / (2 * massH * Navo);
}

__host__ double Cp_from_chi_H(double chi_H, double temperature) {
    double cpH2 = heat_capacity_H2(temperature);
    double cpH  = 20.78603 / (massH * Navo);
    return cpH * chi_H + cpH2 * (1 - chi_H);
}


__device__ double
linear_interpolation_device(double x1, double x2, double y1, double y2, double x) {
    // calculates value of y at x on interval (x1, x2), given y1(x1), y2(x2)
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

__device__ int locate_min_i_device(double *array_d, int N, double val) {

    int id = -1;
    if (val >= array_d[N - 1]) {
        id = N - 2;
    }
    else if (val < array_d[0]) {
        id = 0;
    }
    else {
        for (int j = 1; j < N; j++) {
            if (val >= array_d[j - 1] && val < array_d[j]) {
                id = j - 1;
                break;
            }
        }
    }
    return id;
}

__device__ int locate_max_i_device(double *array_d, int N, double val) {

    int id = -1;
    if (val >= array_d[N - 1]) {
        id = N - 1;
    }
    else if (val < array_d[0]) {
        id = 1;
    }
    else {
        for (int j = 1; j < N; j++) {
            if (val >= array_d[j - 1] && val < array_d[j]) {
                id = j;
                break;
            }
        }
    }
    return id;
}

__device__ double chi_H_equilibrium_device(double *GibbsT,
                                           double *GibbsdG,
                                           int     GibbsN,
                                           double  temperature,
                                           double  pressure) {
    int imin, imax;
    imin         = locate_min_i_device(GibbsT, GibbsN, temperature);
    imax         = locate_max_i_device(GibbsT, GibbsN, temperature);
    double dGnew = linear_interpolation_device(
        GibbsT[imin], GibbsT[imax], GibbsdG[imin], GibbsdG[imax], temperature);
    double kprime = pressure / 1e5 * exp(2000 * dGnew / RUNIV / temperature);
    double n_H    = (-1.0 + sqrt(1 + 8 * kprime)) / (4 * kprime);
    double X_H    = 2 * n_H / (n_H + 1);
    double mmean  = 2 * massH - massH * X_H;
    return massH / mmean
           * X_H; //hmm, in this simple case, chi_H = n_H... wonder if there's any reason to generalize
}

__device__ double Rd_from_chi_H_device(double chi_H) {
    return kBoltz * chi_H / massH + kBoltz * (1 - chi_H) / (2 * massH);
}

__device__ double heat_capacity_H2_device(double temperature) {
    // heat capacity at constant pressure of molecular hydrogen from chase 1998 (nist database)
    double A, B, C, D, E, t = temperature / 1000.0;

    if (temperature < 1000) {
        A = 33.066178;
        B = -11.363417;
        C = 11.432816;
        D = -2.772874;
        E = -0.158558;
    }
    else if (temperature >= 1000 && temperature < 2500) {
        A = 18.563083;
        B = 12.257357;
        C = -2.859786;
        D = 0.268238;
        E = 1.977990;
    }
    else {
        A = 43.423560;
        B = -4.293079;
        C = 1.272428;
        D = -0.096876;
        E = -20.533862;
    }
    return (A + B * t + C * pow(t, 2) + D * pow(t, 3) + E * pow(t, -2)) / (2 * massH * Navo);
}

__device__ double Cp_from_chi_H_device(double chi_H, double temperature) {
    double cpH2 = heat_capacity_H2_device(temperature);
    double cpH  = 20.78603 / (massH * Navo);
    return cpH * chi_H + cpH2 * (1 - chi_H);
}

__global__ void update_temperature_Rd_Cp(double *temperature_d,
                                         double *Rd_d,
                                         double *Cp_d,
                                         double *pressure_d,
                                         double *Rho_d,
                                         double *GibbsT_d,
                                         double *GibbsdG_d,
                                         int     GibbsN,
                                         int     num,
                                         bool *  error) {
    //calculates temperature and Rd(T) from P and Rho using Ridders' method
    //find zeros of function f(T) = P/rho - R*T
    //also update Cp(T)

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double T1 = 0, T2 = 6000, T3, T4 = 6000, T4tmp = 0;
        double f1, f2, f3, f4; // store f(T) = P/rho - R*T
        double eps = 1e-8;
        double p   = pressure_d[id * nv + lev];
        double rho = Rho_d[id * nv + lev];
        double chi_H_tmp;

        //get atomic H mixing ratio at T1 and T2, and evaluate f1, f2
        chi_H_tmp = chi_H_equilibrium(GibbsT_d, GibbsdG_d, GibbsN, T1, p);
        f1        = p / rho - Rd_from_chi_H(chi_H_tmp) * T1;
        chi_H_tmp = chi_H_equilibrium(GibbsT_d, GibbsdG_d, GibbsN, T2, p);
        f2        = p / rho - Rd_from_chi_H(chi_H_tmp) * T2;
        if (f1 * f2 < 0) {
            while (fabs(T4 - T4tmp) > eps) {
                T4tmp = T4;
                T3    = 0.5 * (T1 + T2); //bisect
                // get atomic H mixing ratio at T3 and evaluate f3
                chi_H_tmp = chi_H_equilibrium(GibbsT_d, GibbsdG_d, GibbsN, T3, p);
                f3        = p / rho - Rd_from_chi_H(chi_H_tmp) * T3;
                // false position step for T4
                T4 = T3 + (T3 - T1) * copysignf(1, f1) * f3 / sqrt(pow(f3, 2) - f1 * f2);
                // get atomic H mixing ratio at T4 and evaluate f4
                chi_H_tmp = chi_H_equilibrium(GibbsT_d, GibbsdG_d, GibbsN, T4, p);
                f4        = p / rho - Rd_from_chi_H(chi_H_tmp) * T4;
                if (f4 * f1 > 0) { //T4 is left of the root
                    T1 = T4;
                    f1 = f4;
                }
                else { //T4 is right of the root
                    T2 = T4;
                    f2 = f4;
                }
            }
        }
        else {
            printf("f1=%f,f2=%f", f1, f2);
            *error = true;
        }
        temperature_d[id * nv + lev] = T4;
        Rd_d[id * nv + lev]          = Rd_from_chi_H(chi_H_tmp);
        Cp_d[id * nv + lev]          = Cp_from_chi_H(chi_H_tmp, T4);
        // if (id == 0) {
        //     printf("lev = %d, p = %f, rho = %f, T = %f, Rd = %f, Cp = %f\n",
        //            lev,
        //            p,
        //            rho,
        //            T4,
        //            Rd_d[id * nv + lev],
        //            Cp_d[id * nv + lev]);
        // }
    }
}
