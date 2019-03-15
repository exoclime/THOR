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
// Description: calculates energy, mass, and momentum as a function of time
//
// Method:
//
// Known limitations: None.
//
// Known issues: None.
//
// If you use this code please cite the following references:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
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
#pragma once

__global__ void CalcTotEnergy(double *Etotal_d,
                              double *GlobalE_d,
                              double *Mh_d,
                              double *W_d,
                              double *Rho_d,
                              double *temperature_d,
                              double  Gravit,
                              double  Cp,
                              double  Rd,
                              double  A,
                              double *Altitude_d,
                              double *Altitudeh_d,
                              double *lonlat_d,
                              double *areasT,
                              double *func_r_d,
                              int     num,
                              bool    DeepModel);

__global__ void CalcAngMom(double *AngMomx_d,
                           double *AngMomy_d,
                           double *AngMomz_d,
                           double *GlobalAMx_d,
                           double *GlobalAMy_d,
                           double *GlobalAMz_d,
                           double *Mh_d,
                           double *Rho_d,
                           double  A,
                           double  Omega,
                           double *Altitude_d,
                           double *Altitudeh_d,
                           double *lonlat_d,
                           double *areasT,
                           double *func_r_d,
                           int     num,
                           bool    DeepModel);

__global__ void CalcMass(double *Mass_d,
                         double *GlobalMass_d,
                         double *Rho_d,
                         double  A,
                         double *Altitudeh_d,
                         double *lonlat_d,
                         double *areasT,
                         int     num,
                         bool    DeepModel);

__global__ void CalcEntropy(double *Entropy_d,
                            double *pressure_d,
                            double *temperature_d,
                            double  Cp,
                            double  Rd,
                            double  A,
                            double  P_Ref,
                            double *Altitude_d,
                            double *Altitudeh_d,
                            double *lonlat_d,
                            double *areasT,
                            double *func_r_d,
                            int     num,
                            bool    DeepModel);
