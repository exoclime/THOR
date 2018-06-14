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

// Diffusion
#define HyDiff_default       true    // Hyper-diffusion
#define DivDampP_default     true    // Divergence-damping

// Model options
#define NonHydro_default  true       // Non-hydrostatic parameter
#define DeepModel_default true       // Deep atmosphere

// Initial conditions
#define rest_default   true          // Starting from rest

// Benchmark test
#define hstest_default 1          // Held-Suarez test for Earth == 1
//  HS test for shallow hot Jupiter == 3
//  HS test for tidally locked Earth == 2
//  No HS test == 0

// GPU ID
#define GPU_ID_N_default 0        // Set GPU ID number

// Output
#define n_out_default 1000         // Print output every n_out steps
