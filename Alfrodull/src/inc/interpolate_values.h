// ==============================================================================
// This file is part of Alfrodull.
//
//     Alfrodull is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     Alfrodull is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     Alfrodull directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
// Method: Helios Two Stream algorithm
//
//
// Known limitations: - Runs in a single GPU.
//
// Known issues: None
//
//
// Code contributors: Urs Schroffenegger, Matej Malik
//
// History:
// Version Date       Comment
// ======= ====       =======
// 1.0     2020-07-15 First version
//
//
////////////////////////////////////////////////////////////////////////

__global__ void interpolate_temperature(double* tlay, double* tint, int numinterfaces);

// interpolates the Planck function for the layer temperatures from the pre-tabulated values
__global__ void planck_interpol_layer(double* temp,
                                      double* planckband_lay,
                                      double* planck_grid,
                                      double* starflux,
                                      bool    null_planck_function,
                                      bool    realstar,
                                      int     numlayers,
                                      int     nwave,
                                      int     dim,
                                      int     step);


// interpolates the Planck function for the interface temperatures from the pre-tabulated values
__global__ void planck_interpol_interface(double* temp,
                                          double* planckband_int,
                                          double* planck_grid,
                                          bool    null_planck_function,
                                          int     numinterfaces,
                                          int     nwave,
                                          int     dim,
                                          int     step);

// interpolate layer and interface opacities from opacity table
__global__ void interpolate_opacities(double* temp,
                                      double* opactemp,
                                      double* press,
                                      double* opacpress,
                                      double* ktable,
                                      double* opac,
                                      double* crosstable,
                                      double* scat_cross,
                                      int     npress,
                                      int     ntemp,
                                      int     ny,
                                      int     nbin,
                                      double  opaclimit,
                                      int     temp_num_per_col,
                                      int     press_num_per_col,
                                      int     nlay_or_nint);


// interpolate the mean molecular mass for each layer
__global__ void meanmolmass_interpol(double* temp,
                                     double* opactemp,
                                     double* meanmolmass,
                                     double* opac_meanmass,
                                     double* press,
                                     double* opacpress,
                                     int     npress,
                                     int     ntemp,
                                     int     temp_num_per_col,
                                     int     press_num_per_col,
                                     int     ninterface);


// interpolate kappa for each layer
__global__ void kappa_interpol(double* temp,
                               double* entr_temp,
                               double* press,
                               double* entr_press,
                               double* kappa,
                               double* opac_kappa,
                               int     entr_npress,
                               int     entr_ntemp,
                               int     nlay_or_nint,
                               double  kappa_kernel_value);
