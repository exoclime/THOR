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
// ESP -  Exoclimes Simulation Platform. (version 1.0)
//
//   
//
// Method: [1] - Builds the icosahedral grid (Mendonca et al., 2016).
//         Iteration:
//         [2] - Calls the dynamical core THOR (Mendonca et al., 2016).
//         [3] - Calls the physical core ProfX.
//         END
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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

using namespace std;

#include "headers/define.h"
#include "headers/grid.h"
#include "headers/planet.h"
#include "headers/esp.h"

#include "debug.h"

#include "config_file.h"

int main (){

    printf("\n Starting ESP!");
    printf("\n\n version 1.0\n");
    printf(" Compiled on " __DATE__ " at " __TIME__ "\n");
    
    printf(" build level: " BUILD_LEVEL "\n");

    // Initial conditions
    config_file config_reader;

    // Set config variables
    
    //  planet // initialises to default earth
    XPlanet Planet;
    
    int timestep = 1000;
    int nsmax = 48000;
    string simulation_ID = "Earth";
    

    config_reader.append_config_var("timestep", timestep, timestep_default);
    config_reader.append_config_var("num_steps", nsmax, nsmax_default);
    config_reader.append_config_var("simulation_ID", simulation_ID, string(Planet.simulation_ID));
    config_reader.append_config_var("radius", Planet.A, Planet.A);
    config_reader.append_config_var("rotation_rate", Planet.Omega, Planet.Omega);
    config_reader.append_config_var("gravitation", Planet.Gravit, Planet.Gravit);
    config_reader.append_config_var("Mmol", Planet.Mmol, Planet.Mmol);
    config_reader.append_config_var("Rd", Planet.Rd, Planet.Rd);
    config_reader.append_config_var("Cp", Planet.Cp, Planet.Cp);
    config_reader.append_config_var("Tmean", Planet.Tmean, Planet.Tmean);
    config_reader.append_config_var("P_ref", Planet.P_Ref, Planet.P_Ref);
    config_reader.append_config_var("Top_altitude", Planet.Top_altitude, Planet.Top_altitude);
    config_reader.append_config_var("Diffc", Planet.Diffc, Planet.Diffc);
    // Read config file
    // TODO: should check some default path
    // TODO: If nothing given, should fall back to 
    // TODO: should be overridable by command line options
    string initial_filename = "ifile/earth.thr";

    if (config_reader.parse_file(initial_filename))
        printf(" Config file %s read\n", initial_filename.c_str());
    else
    {
        printf(" Config file %s reading failed, aborting\n", initial_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    // Test config variables for coherence
    bool config_OK = true;    
    config_OK &= check_greater( "timestep", timestep, 0);
    config_OK &= check_greater( "nsmax", nsmax, 0);
    config_OK &= check_greater( "gravitation", Planet.Gravit, 0.0);
    config_OK &= check_greater( "Mmol", Planet.Mmol, 0.0);
    config_OK &= check_greater( "Rd", Planet.Rd, 0.0);
    config_OK &= check_greater( "T_mean", Planet.Tmean, 0.0);
    config_OK &= check_greater( "P_Ref", Planet.P_Ref, 0.0);
    config_OK &= check_greater( "Top_altitude", Planet.Top_altitude, 0.0);

    if (simulation_ID.length() < 160)
    {
        sprintf(Planet.simulation_ID, "%s", simulation_ID.c_str());
    }
    else
    {
        printf("Bad value for config variable simulation_ID: [%s]\n", simulation_ID.c_str());
        
        config_OK = false;
    }
    
    
    
    if (!config_OK)
    {
        printf("Error in configuration file\n");
        exit(0);
    }
    
//  Set the GPU device.    
    cudaError_t error;
    cudaSetDevice(GPU_ID_N);
    
//
//  Make the icosahedral grid
    Icogrid Grid(sprd               , // Spring dynamics option
                 spring_beta        , // Parameter beta for spring dynamics 
                 glevel             , // Horizontal resolution level
                 vlevel             , // Number of vertical layers
                 Planet.A           , // Planet radius
                 Planet.Top_altitude);// Top of the model's domain

//
//  Define object X.
    ESP X( Grid.point_local   , // First neighbours
           Grid.maps          , // Grid domains
           Grid.lonlat        , // Longitude and latitude of the grid points
           Grid.Altitude      , // Altitudes
           Grid.Altitudeh     , // Altitude at the interfaces between layers
           Grid.nvecoa        , // Normal vectors for diffusion 1
           Grid.nvecti        , // Normal vectors for diffusion 2
           Grid.nvecte        , // Normal vectors for diffusion 3
           Grid.areasT        , // Areas of the main cells
           Grid.areasTr       , // Areas of the triangles
           Grid.div           , // Divergence operator
           Grid.grad          , // Gradient operator
           Grid.func_r        , // Normalised vector
           Grid.nl_region     , // Number of points in one side of a rhombus
           Grid.nr            , // Number of rhombi
           Grid.nv            , // Number of vertical layers
           Grid.nvi           , // Number of interfaces between layers
           Grid.point_num     );// Number of grid points

    USE_BENCHMARK();
    
    BENCH_POINT_GRID(Grid);
    
   printf(" Setting the initial conditions.\n\n");
// Initial conditions
   X.InitialValues(rest         , // Option to start the atmosphere from rest
                   glevel       , // Horizontal resolution level
                   timestep     , // Time-step [s]
                   Planet.A     , // Planet radius [m]
                   Planet.Cp    , // Specific heat capacity [J /(kg K)]
                   Planet.P_Ref , // Reference pressure [Pa]
                   Planet.Gravit, // Gravity [m/s^2]
                   Planet.Omega , // Rotation rate [1/s]
                   Planet.Diffc , // Strength of diffusion 
                   kb           , // Boltzmann constant [J/kg]
                   Planet.Tmean , // Isothermal atmosphere (at temperature Tmean)
                   Planet.Mmol  , // Mean molecular mass of dry air [kg]
                   mu           , // Atomic mass unit [kg]
                   Planet.Rd    );// Gas constant [J/kg/K]

    long startTime = clock();

//
//  PRINTS
//  Device Information
    int ndevices;
    cudaGetDeviceCount(&ndevices);
    for (int i = 0; i < ndevices; ++i) {
        // Get device properties
        printf("\n CUDA Device #%d\n", i);
        cudaDeviceProp devPp;
        cudaGetDeviceProperties(&devPp, i);
        printf(" Name: %s\n",  devPp.name);
        printf("   Total global memory:           %lu\n",  devPp.totalGlobalMem);
        printf("   Total shared memory per block: %lu\n",  devPp.sharedMemPerBlock);
        printf("   Total registers per block:     %d\n" ,  devPp.regsPerBlock);
        printf("   Warp size:                     %d\n" ,  devPp.warpSize);
        printf("   Maximum memory pitch:          %lu\n",  devPp.memPitch);
        printf("   Maximum threads per block:     %d\n" ,  devPp.maxThreadsPerBlock);
        printf("   Clock rate:                    %d\n" ,  devPp.clockRate);
        printf("   Total constant memory:         %lu\n",  devPp.totalConstMem);
        printf("   Number of multiprocessors:     %d\n" ,  devPp.multiProcessorCount);
    }
//
//  Planet conditions
    printf("\n");
    printf(" Planet: %s\n",Planet.simulation_ID);
    printf("   Radius = %f m\n"       , Planet.A     );
    printf("   Omega  = %f s-1\n"     , Planet.Omega );
    printf("   Gravit = %f m/s2\n"    , Planet.Gravit);
    printf("   Mmol   = %f kg\n"      , Planet.Mmol  );
    printf("   Rd     = %f J/(Kg K)\n", Planet.Rd    );
    printf("   Cp     = %f J/(Kg K)\n", Planet.Cp    );
    printf("   Tmean  = %f K\n"       , Planet.Tmean );
//
//  Numerical Methods
    printf("\n");
    printf(" THOR\n");
    printf("   ********** \n");
    printf("   Grid - Icosahedral\n");
    printf("   Glevel          = %d.\n", glevel);
    printf("   Spring dynamics = %d.\n",sprd);
    printf("   Beta            = %f.\n", spring_beta);
    printf("   Resolution      = %f deg.\n", (180/M_PI)*sqrt(2*M_PI/5)/pow(2,glevel));
    printf("   Vertical layers = %d.\n",Grid.nv);
    printf("   ********** \n");
    printf("   Split-Explicit / HE-VI \n");
    printf("   FV = Central finite volume \n");
    printf("   Time integration =  %d s.\n", nsmax*timestep);
    printf("   Large time-step  =  %d s.\n", timestep);
    printf("    \n");
    
//
//  Writes initial conditions
    double simulation_time = 0.0;
    X.Output(0                   ,
             Planet.Cp           , // Specific heat capacity [J/(Kg K)]
             Planet.Rd           , // Gas constant [J/(Kg K)]
             Planet.Omega        , // Rotation rate [s-1]
             Planet.Gravit       , // Gravitational acceleration [m/s2]
             Planet.Mmol         , // Mean molecular mass of dry air [kg]
             Planet.P_Ref        , // Reference surface pressure [Pa]
             Planet.Top_altitude , // Top of the model's domain [m]
             Planet.A            , // Planet Radius [m]
             Planet.simulation_ID, // Simulation ID (e.g., "Earth")
             simulation_time     );// Time of the simulation [s]

//
//  Starting model Integration.
    printf(" Starting the model integration.\n\n");

//
//  Main loop. nstep is the current step of the integration and nsmax the maximum
//  number of steps in the integration.
    for(int nstep = 1; nstep <= nsmax; ++nstep){
        X.current_step = nstep;
        
//
//      Dynamical Core Integration (THOR)
        X.Thor (timestep     , // Time-step [s]
                HyDiff       , // Hyperdiffusion option
                DivDampP     , // Divergence-damping option
                Planet.Omega , // Rotation rate [1/s]
                Planet.Cp    , // Specific heat capacity [J/kg/K]
                Planet.Rd    , // Gas constant [J/kg/K]
                Planet.Mmol  , // Mean molecular mass of dry air [kg]
                mu           , // Atomic mass unit [kg]
                kb           , // Boltzmann constant [J/K]
                Planet.P_Ref , // Reference pressure [Pa]
                Planet.Gravit, // Gravity [m/s^2]
                Planet.A     , // Planet radius [m]
                NonHydro     , // Non-hydrostatic option
                DeepModel    );// Deep model option

//
//     Physical Core Integration (ProfX)
       X.ProfX(planetnumber , // Planet ID 
               nstep        , // Step number
               hstest       , // Held-Suarez test option
               timestep     , // Time-step [s]
               Planet.Omega , // Rotation rate [1/s]
               Planet.Cp    , // Specific heat capacity [J/kg/K]
               Planet.Rd    , // Gas constant [J/kg/K]
               Planet.Mmol  , // Mean molecular mass of dry air [kg]
               mu           , // Atomic mass unit [kg]
               kb           , // Boltzmann constant [J/K]
               Planet.P_Ref , // Reference pressure [Pa]
               Planet.Gravit, // Gravity [m/s^2]
               Planet.A     );// Planet radius [m]

//
//      Prints output every nout steps
        if(nstep % n_out == 0) {
            X.CopyToHost();
            X.Output(nstep/n_out         ,
                     Planet.Cp           , // Specific heat capacity [J/(Kg K)]
                     Planet.Rd           , // Gas constant [J/(Kg K)]
                     Planet.Omega        , // Rotation rate [s-1]
                     Planet.Gravit       , // Gravitational acceleration [m/s2]
                     Planet.Mmol         , // Mean molecular mass of dry air [kg]
                     Planet.P_Ref        , // Reference surface pressure [Pa] 
                     Planet.Top_altitude , // Top of the model's domain [m]
                     Planet.A            , // Planet radius [m]
                     Planet.simulation_ID, // Planet ID
                     simulation_time     );// Simulation time [s]
        }
        printf("\n Time step number = %d || Time = %f days.", nstep, nstep*timestep/86400.);
    }
//
//  Prints the duration of the integration.
    long finishTime = clock();
    cout << "\n\n Integration time = " << (finishTime - startTime)/CLOCKS_PER_SEC 
                                       << " seconds" << endl;

//
//  Checks for errors in the device.
    cudaDeviceSynchronize();
    error = cudaGetLastError();
    printf("\n\n CudaMalloc error = %d = %s",error, cudaGetErrorString(error));
//
//  END OF THE ESP
    printf("\n\n End of the ESP!\n\n");

    return 0;
}
