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
// Description: dry convective adjustment scheme
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

#include <math.h>

#define Gibbs_Filename "src/headers/phy/GibbsH.txt"
#define RUNIV 8.3144621
#define massH 1.67372e-27
#define kBoltz 1.38064852e-23
#define Navo 6.0221409e23

__host__ bool ESP::read_in_gibbs_H() {
    FILE *infile;
    int   nlines = 61; //lines in GibbsH file

    infile  = fopen(Gibbs_Filename, "r");
    GibbsT  = (double *)malloc(nlines * sizeof(double));
    GibbsdG = (double *)malloc(nlines * sizeof(double));

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

__host__ double linear_interpolation(double x1, double x2, double y1, double y2, double x) {
    // calculates value of y at x on interval (x1, x2), given y1(x1), y2(x2)
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

__host__ int locate_min_i(double *array_d, int N, double val) {

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

__host__ int locate_max_i(double *array_d, int N, double val) {

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

__host__ double kprime(double *T, double *dG, double temperature, double pressure) {
    //calculates kprime from Equation 48 in Heng, Lyons, and Tsai 2016
    int imin, imax;
    imin         = locate_min_i(T, 61, temperature);
    imax         = locate_max_i(T, 61, temperature);
    double dGnew = linear_interpolation(T[imin], T[imax], dG[imin], dG[imax], temperature);
    return pressure / 1e5 * exp(2000 * dGnew / RUNIV / temperature);
}

__host__ double ESP::chi_H_equilibrium(double temperature, double pressure) {
    int imin, imax;
    imin = locate_min_i(GibbsT, 61, temperature);
    imax = locate_max_i(GibbsT, 61, temperature);
    double dGnew =
        linear_interpolation(GibbsT[imin], GibbsT[imax], GibbsdG[imin], GibbsdG[imax], temperature);
    double kprime = pressure / 1e5 * exp(2000 * dGnew / RUNIV / temperature);
    double n_H    = (-1.0 + sqrt(1 + 8 * kprime)) / (4 * kprime);
    double X_H    = 2 * n_H / (n_H + 1);
    double mmean  = 2 * massH - massH * X_H;
    return massH / mmean
           * X_H; //hmm, in this simple case, chi_H = n_H... wonder if there's any reason to generalize
}

__host__ double Rd_from_chi_H(double chi_H) {
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
