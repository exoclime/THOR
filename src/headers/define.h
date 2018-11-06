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
// Description: Defines the main model's parameters
//
// Method: -
//
// Known limitations: None
//
//
// Known issues: None
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

// Integration time
#define nsmax_default        48000 // Number of time steps
#define timestep_default     1800  // Time step  [seconds]

// Grid
#define sprd_default   true       // Spring dynamics2
#define glevel_default 4          // Horizontal resolution level.
#define spring_beta_default 1.15  // Parameter beta for spring dynamics
#define vlevel_default 32         // Number of vertical layers

// Sponge layer
#define nlat_default   20         // Number of latitude rings for mean zonal wind (sponge layer)
#define Rv_sponge_default 1e-4    // Maximum damping (top of model)
#define ns_sponge_default 0.75    // Lowest level of sponge layer
#define shrink_sponge_default false       // shrink sponge layer after some time
#define t_shrink_default 500.0            // shrink sponge after this many days

// Diffusion
#define HyDiff_default       true    // Hyper-diffusion
#define DivDampP_default     true    // Divergence-damping

// Model options
#define NonHydro_default  true       // Non-hydrostatic parameter
#define DeepModel_default true       // Deep atmosphere
#define SpongeLayer_default false    // use sponge layer at top of model

// Initial conditions
#define rest_default   true          // Starting from rest
#define initial_conditions_default   "ifile/esp_initial.h5" // start
                                                            // from
                                                            // this
                                                            // initial
                                                            // conditions
                                                            // file

// Benchmark test
#define core_benchmark_default 1          // Held-Suarez test for Earth == 1
//  HS test for shallow hot Jupiter == 3
//  HS test for tidally locked Earth == 2
//  No HS test == 0

// GPU ID
#define GPU_ID_N_default 0        // Set GPU ID number

// Output
#define n_out_default 1000         // Print output every n_out steps


#define output_path_default  "results"     // Output directory

#define gcm_off_default false      //turns off fluid dynamical core for debugging physics
#define TPprof_default 0           //initial TP profile (0 = isothermal, 1 = guillot)

#define conservation_default false //output energy, mass, angular momentum, etc

#define vulcan_default 0   //use chemical kinetics scheme
