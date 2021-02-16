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

// Integration time
#define nsmax_default 48000   // Number of time steps
#define timestep_default 1800 // Time step  [seconds]

// Grid
#define sprd_default true        // Spring dynamics2
#define glevel_default 4         // Horizontal resolution level.
#define spring_beta_default 1.15 // Parameter beta for spring dynamics
#define vlevel_default 32        // Number of vertical layers

// Sponge layer
#define RayleighSponge_default false  // use sponge layer at top of model (rayleigh drag)
#define RayleighSpongeT_default false // use T term in sponge layer (rayleigh drag)
#define nlat_bins_default 20          // Number of latitude rings for mean zonal wind (sponge layer)
#define Ruv_sponge_default 1e-4       // Maximum rayleigh damping (top of model)
#define Rw_sponge_default 1e-4        // Maximum tayleigh damping (top of model)
#define RT_sponge_default 1e-4        // Maximum damping (top of model)
#define ns_ray_sponge_default 0.75    // Lowest level of sponge layer
#define raysp_calc_mode_default "imp" // method for updating winds from rayleigh sponge
#define damp_uv_to_mean_default true  // damps uv toward zonal mean (otherwise to zero)
#define damp_w_to_mean_default true   // damps w toward zonal mean (otherwise to zero)

#define DiffSponge_default false    // use diffusive sponge
#define ns_diff_sponge_default 0.75 // Lowest level of sponge layer
#define Dv_sponge_default 0.01      // Maximum diff damping (top of model)
#define order_diff_sponge_default 2 // order of diffusive sponge

#define shrink_sponge_default false // shrink sponge layer after some time
#define t_shrink_default 144000     // shrink sponge after this many time steps

// Diffusion
#define HyDiff_default true // Hyper-diffusion
#define HyDiffOrder_default 4
#define DivDampP_default true    // Divergence-damping
#define VertHyDiff_default false // Vertical Hyper-diffusion

// Model options
#define NonHydro_default true    // Non-hydrostatic parameter
#define DeepModel_default true   // Deep atmosphere
#define output_mean_default true // output mean quantities

// Initial conditions
#define rest_default true                                 // Starting from rest
#define initial_conditions_default "ifile/esp_initial.h5" // start from this initial conditions file

// Benchmark test
#define core_benchmark_default "HeldSuarez" // Held-Suarez test for Earth == "HeldSuarez"

// GPU ID
#define GPU_ID_N_default 0 // Set GPU ID number

// Output
#define n_out_default 1000                // Print output every n_out steps
#define custom_global_n_out_default false // Print global diagnostics every n_out or global_n_out
#define global_n_out_default 1000         // Print global diagnostics every global_n_out steps
#define custom_log_n_out_default false    // Print log info every n_out steps or log_n_out
#define log_n_out_default 1000            // Print log info every log_n_out steps

#define output_path_default "results" // Output directory

#define gcm_off_default false //turns off fluid dynamical core for debugging physics

#define globdiag_default false //output energy, mass, angular momentum, etc

#define conv_adj_default false // use convective adjustment scheme
// number of times to execute per time step
// (repeats entire algorithm if > 1)
#define conv_adj_iter_default 1

#define init_PT_profile_default "isothermal"
#define kappa_lw_default 0.002 // m^2 kg^-1
#define kappa_sw_default 0.001 // m^2 kg^-1
#define Tint_default 100.0
#define f_lw_default 0.5
#define bv_freq_default 0.01

#define uh_thermo_default "none"

#define uh_heating_default "none"

#define thermo_equation_default "entropy"

#define vert_refined_default false
#define n_bl_layers_default 9

enum benchmark_types {
    NO_BENCHMARK         = 0,
    HELD_SUAREZ          = 1,
    TIDALLY_LOCKED_EARTH = 2,
    SHALLOW_HOT_JUPITER  = 3,
    DEEP_HOT_JUPITER     = 4,
    JET_STEADY           = 5,
    ACOUSTIC_TEST        = 6,
    GWAVE_TEST           = 7
};

enum init_PT_profile_types { ISOTHERMAL = 0, GUILLOT = 1, CONSTBV = 2 };

enum uh_thermo_types { NO_UH_THERMO = 0, VARY_R_CP = 1, FULL = 2 };

enum uh_heating_types { NO_UH_HEATING = 0, QUASI_EQL = 1, RELAX_CHEM = 2 };

enum raysp_calc_mode_types { IMP = 0, EXP1 = 1, EXP3 = 2 };

enum thermo_equation_types { ENTROPY = 0, ENERGY = 1 }; //might add pressure?
