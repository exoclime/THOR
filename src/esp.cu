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
// ESP -  Exoclimes Simulation Platform. (version 2.0)
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

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <iomanip>
#include <sstream>

#include "esp.h"
#include "headers/define.h"
#include "headers/grid.h"
#include "headers/planet.h"

#include "cmdargs.h"
#include "config_file.h"
#include "directories.h"
#include "iteration_timer.h"

#include "binary_test.h"

#include "phy_modules.h"

#include "debug.h"
#ifdef BENCHMARKING
#    warning "Compiling with benchmarktest enabled"
#endif

#include <csignal>

using std::string;


enum e_sig {
    ESIG_NOSIG   = 0,
    ESIG_SIGTERM = 1,
    ESIG_SIGINT  = 2
};


volatile sig_atomic_t caught_signal = ESIG_NOSIG;


void sigterm_handler(int sig) {
    caught_signal = ESIG_SIGTERM;
    printf("SIGTERM caught, trying to exit gracefully.\n");
}

void sigint_handler(int sig) {
    caught_signal = ESIG_SIGINT;
    printf("SIGINT caught, trying to exit gracefully\n");
}


std::string duration_to_str(double delta) {
    unsigned int days = delta / (24 * 3600);
    delta -= days * 24.0 * 3600.0;

    unsigned int hours = delta / 3600;
    delta -= hours * 3600;

    unsigned int minutes = delta / 60;
    delta -= minutes * 60;

    unsigned int       seconds = delta;
    std::ostringstream str;

    if (days != 0)
        str << days << "d ";
    if (hours != 0)
        str << hours << "h ";
    if (minutes != 0)
        str << minutes << "m ";
    str << seconds << "s";

    return str.str();
}

void get_cuda_mem_usage(size_t& total_bytes, size_t& free_bytes) {
    // show memory usage of GPU
    cudaError_t cuda_status = cudaMemGetInfo(&free_bytes, &total_bytes);

    if (cudaSuccess != cuda_status) {
        printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
    }
}


int main(int argc, char** argv) {
    //*****************************************************************
    // Parse input arguments
    cmdargs argparser("esp", "run THOR GCM simulation");

    // define arguments
    // format:
    // * short form: e.g. "i" for -i
    // * long form: e.g. "int" for --int
    // * default value, defines type of variable: e.g. -1 for int,
    //   value -1
    // * help string
    // argparser.add_arg("i", "int", -1, "This is an int");
    // argparser.add_arg("b", "bool", false, "this is a bool");
    // argparser.add_arg("d", "double", 1e-5, "this is a double");
    // argparser.add_arg("s", "string", string("value"), "this is a string");

    string default_config_filename("ifile/earth.thr");

    argparser.add_positional_arg(default_config_filename, "config filename");
    argparser.add_arg("g", "gpu_id", 0, "GPU_ID to run on");
    argparser.add_arg("o", "output_dir", string("results"), "results directory to store output");

    argparser.add_arg("i", "initial", string("initialfilename"), "start from this initial condition instead of rest");
    argparser.add_arg("c", "continue", string("continuefilename"), "continue simulation from this output file");
    argparser.add_arg("N", "numsteps", 48000, "number of steps to run");
    argparser.add_arg("w", "overwrite", true, "Force overwrite of output file if they exist");

    argparser.add_arg("b", "batch", true, "Run as batch");


    // Parse arguments, exit on fail
    bool parse_ok = argparser.parse(argc, argv);
    if (!parse_ok) {
        printf("error parsing arguments\n");

        argparser.print_help();
        exit(-1);
    }


    // get config file path if present in argument
    string config_filename = "";
    argparser.get_positional_arg(config_filename);

    //*****************************************************************
    printf("\n Starting ESP!");
    printf("\n\n version 1.0\n");
    printf(" Compiled on %s at %s.\n", __DATE__, __TIME__);

    printf(" build level: %s\n", BUILD_LEVEL);


    //*****************************************************************
    // Initial conditions
    config_file config_reader;


    //*****************************************************************
    // Setup config variables

    // Integration time
    int timestep = 1000;
    int nsmax    = 48000;

    config_reader.append_config_var("timestep", timestep, timestep_default);
    config_reader.append_config_var("num_steps", nsmax, nsmax_default);

    // Planet
    // initialises to default earth
    XPlanet Planet;
    string  simulation_ID = "Earth";

    config_reader.append_config_var("simulation_ID", simulation_ID, string(Planet.simulation_ID));
    config_reader.append_config_var("radius", Planet.A, Planet.A);
    config_reader.append_config_var("rotation_rate", Planet.Omega, Planet.Omega);
    config_reader.append_config_var("gravitation", Planet.Gravit, Planet.Gravit);
    config_reader.append_config_var("Rd", Planet.Rd, Planet.Rd);
    config_reader.append_config_var("Cp", Planet.Cp, Planet.Cp);
    config_reader.append_config_var("Tmean", Planet.Tmean, Planet.Tmean);
    config_reader.append_config_var("P_ref", Planet.P_Ref, Planet.P_Ref);
    config_reader.append_config_var("Top_altitude", Planet.Top_altitude, Planet.Top_altitude);
    config_reader.append_config_var("Diffc", Planet.Diffc, Planet.Diffc);

    // grid
    bool   spring_dynamics = true;
    int    glevel          = 4;
    double spring_beta     = 1.15;
    int    vlevel          = 32;

    config_reader.append_config_var("spring_dynamics", spring_dynamics, sprd_default);
    config_reader.append_config_var("glevel", glevel, glevel_default);
    config_reader.append_config_var("spring_beta", spring_beta, spring_beta_default);
    config_reader.append_config_var("vlevel", vlevel, vlevel_default);

    // Diffusion
    bool HyDiff   = true;
    bool DivDampP = true;
    config_reader.append_config_var("HyDiff", HyDiff, HyDiff_default);
    config_reader.append_config_var("DivDampP", DivDampP, DivDampP_default);

    // Model options
    bool   NonHydro      = true;
    bool   DeepModel     = true;
    bool   SpongeLayer   = false;
    int    nlat          = 20;
    double Rv_sponge     = 1e-4;
    double ns_sponge     = 0.75;
    bool   shrink_sponge = false;
    double t_shrink      = 500;

    config_reader.append_config_var("NonHydro", NonHydro, NonHydro_default);
    config_reader.append_config_var("DeepModel", DeepModel, DeepModel_default);
    config_reader.append_config_var("SpongeLayer", SpongeLayer, SpongeLayer_default);
    config_reader.append_config_var("nlat", nlat, nlat_default);
    config_reader.append_config_var("Rv_sponge", Rv_sponge, Rv_sponge_default);
    config_reader.append_config_var("ns_sponge", ns_sponge, ns_sponge_default);
    config_reader.append_config_var("shrink_sponge", shrink_sponge, shrink_sponge_default);
    config_reader.append_config_var("t_shrink", t_shrink, t_shrink_default);

    // Initial conditions
    // rest supersedes initial condition entry,
    // but if initial condition set from  command line, it overrides
    // rest variable
    bool rest = true;
    config_reader.append_config_var("rest", rest, rest_default);

    string initial_conditions = "initialfilename.h5";
    config_reader.append_config_var("initial", initial_conditions, string(initial_conditions_default));

    // Benchmark test
    string core_benchmark_str("HeldSuarez");
    config_reader.append_config_var("core_benchmark", core_benchmark_str, string(core_benchmark_default));

    int conv_adj = 1;
    config_reader.append_config_var("conv_adj", conv_adj, conv_adj_default);

    int GPU_ID_N = 0;
    config_reader.append_config_var("GPU_ID_N", GPU_ID_N, GPU_ID_N_default);

    int n_out = 1000;
    config_reader.append_config_var("n_out", n_out, n_out_default);

    string output_path = "results";
    config_reader.append_config_var("results_path", output_path, string(output_path_default));

    bool gcm_off = false;
    config_reader.append_config_var("gcm_off", gcm_off, gcm_off_default);

    int TPprof = 0;
    config_reader.append_config_var("TPprof", TPprof, TPprof_default);

    bool conservation = false;
    config_reader.append_config_var("conservation", conservation, conservation_default);
    //*****************************************************************
    // read configs for modules
    phy_modules_generate_config(config_reader);


    //*****************************************************************
    // Read config file
    printf("\n");

    if (config_reader.parse_file(config_filename))
        printf(" Config file %s read\n", config_filename.c_str());
    else {
        printf(" Config file %s reading failed, aborting\n", config_filename.c_str());
        exit(-1);
    }

    //*****************************************************************
    // Override config file variables from command line if set
    int GPU_ID_N_arg;
    if (argparser.get_arg("gpu_id", GPU_ID_N_arg))
        GPU_ID_N = GPU_ID_N_arg;

    string output_dir_arg;

    if (argparser.get_arg("output_dir", output_dir_arg))
        output_path = output_dir_arg;

    string inital_conditions_arg;
    bool   initial_condition_arg_set = false;


    if (argparser.get_arg("initial", inital_conditions_arg)) {
        rest                      = false;
        initial_conditions        = inital_conditions_arg;
        initial_condition_arg_set = true;

        if (!path_exists(initial_conditions)) {
            printf("Initial conditions file \"%s\" not found\n", initial_conditions.c_str());
            exit(-1);
        }
    }

    bool run_as_batch_arg = false;
    bool run_as_batch     = false;
    if (argparser.get_arg("batch", run_as_batch_arg)) {
        run_as_batch = run_as_batch_arg;
    }

    string continue_filename = "";
    bool   continue_sim      = false;


    if (argparser.get_arg("continue", continue_filename)) {
        rest         = false;
        continue_sim = true;

        if (run_as_batch) {
            printf("--continue and --batch options set, options are exclusive, must set only one\n");

            exit(-1);
        }


        if (initial_condition_arg_set) {
            printf("--continue and --initial options set, options are exclusive, must set only one\n");

            exit(-1);
        }

        if (!path_exists(continue_filename)) {
            printf("Continuation start condition file \"%s\" not found\n", continue_filename.c_str());
            exit(-1);
        }

        initial_conditions = continue_filename;
    }

    bool force_overwrite_arg = false;
    bool force_overwrite     = false;

    if (argparser.get_arg("overwrite", force_overwrite_arg))
        force_overwrite = force_overwrite_arg;


    int nsmax_arg;
    if (argparser.get_arg("numsteps", nsmax_arg))
        nsmax = nsmax_arg;


    // Test config variables for coherence
    bool config_OK = true;
    config_OK &= check_greater("timestep", timestep, 0);
    config_OK &= check_greater("nsmax", nsmax, 0);
    config_OK &= check_greater("gravitation", Planet.Gravit, 0.0);
    config_OK &= check_greater("Rd", Planet.Rd, 0.0);
    config_OK &= check_greater("T_mean", Planet.Tmean, 0.0);
    config_OK &= check_greater("P_Ref", Planet.P_Ref, 0.0);
    config_OK &= check_greater("Top_altitude", Planet.Top_altitude, 0.0);

    config_OK &= check_range("glevel", glevel, 3, 8);
    config_OK &= check_greater("vlevel", vlevel, 0);

    config_OK &= check_greater("GPU_ID_N", GPU_ID_N, -1);
    config_OK &= check_greater("n_out", n_out, 0);

    config_OK &= check_range("TPprof", TPprof, -1, 2);

    if (simulation_ID.length() < 160) {
        sprintf(Planet.simulation_ID, "%s", simulation_ID.c_str());
    }
    else {
        printf("Bad value for config variable simulation_ID: [%s]\n", simulation_ID.c_str());

        config_OK = false;
    }

    // Check core benchmark string
    benchmark_types core_benchmark = HELD_SUAREZ;

    if (core_benchmark_str == "NoBenchmark") {
        core_benchmark = NO_BENCHMARK;
        config_OK &= true;
    }
    else if (core_benchmark_str == "HeldSuarez") {
        core_benchmark = HELD_SUAREZ;
        config_OK &= true;
    }
    else if (core_benchmark_str == "ShallowHotJupiter") {
        core_benchmark = SHALLOW_HOT_JUPITER;
        config_OK &= true;
    }
    else if (core_benchmark_str == "DeepHotJupiter") {
        core_benchmark = DEEP_HOT_JUPITER;
        config_OK &= true;
    }
    else if (core_benchmark_str == "TidallyLockedEarth") {
        core_benchmark = TIDALLY_LOCKED_EARTH;
        config_OK &= true;
    }
    else if (core_benchmark_str == "JetSteady") {
        core_benchmark = JET_STEADY;
        config_OK &= true;
    }
    else {
        printf("core_benchmark config item not recognised: [%s]\n", core_benchmark_str.c_str());
        config_OK &= false;
    }


    if (!config_OK) {
        printf("Error in configuration file\n");
        exit(-1);
    }
    //*****************************************************************
#ifdef BENCHMARKING
    string output_path_ref = (path(output_path) / string("ref")).to_string();
#    ifdef BENCH_POINT_COMPARE
    output_path = (path(output_path) / string("compare")).to_string();
#    endif // BENCH_POINT_COMPARE
#    ifdef BENCH_POINT_WRITE
    output_path = (path(output_path) / string("write")).to_string();
#    endif // BENCH_POINT_WRITE
#endif     // BENCHMARKING


    // check output config directory
    if (!create_output_dir(output_path)) {
        printf("Error creating output result directory: %s\n",
               output_path.c_str());
        exit(-1);
    }

    //*****************************************************************
    log_writer logwriter(Planet.simulation_ID, output_path);

    // Batch mode handling
    if (run_as_batch) {
        printf("Starting in batch mode.\n");


        // Get last written file from
        int    last_file_number      = 0;
        int    last_iteration_number = 0;
        string last_file             = "";
        bool   has_last_file         = false;

        try {
            has_last_file = logwriter.check_output_log(last_file_number, last_iteration_number, last_file);
        } catch (const std::exception& e) {
            printf("[%s:%d] error while checking output log: %s.\n", __FILE__, __LINE__, e.what());
            exit(-1);
        }

        if (has_last_file) {
            path o(output_path);
            o /= last_file;

            if (path_exists(o.to_string())) {
                printf("continuing batch run from file %s\n", last_file.c_str());

                // we want to continue the simulation
                rest         = false;
                continue_sim = true;

                // overwrite the results files we'd encounter
                force_overwrite = true;

                // reload the last file we found as initial conditions
                initial_conditions = o.to_string();

                logwriter.open_output_log_for_write(true /*open in append mode */);
                logwriter.prepare_conservation_file(true);
                logwriter.prepare_diagnostics_file(true);
            }
            else {
                printf("Did not find last saved file that should exist.\n");
                exit(-1);
            }
        }
        else {
            printf("No batch file found, initialise simulation.\n");
            // we don't have an output file, start from scratch, reinitialising outputs
            logwriter.open_output_log_for_write(false /*open in non append mode */);
            logwriter.prepare_conservation_file(false);
            logwriter.prepare_diagnostics_file(false);
        }
    }
    else {
        printf("Opening result output file.\n");
        logwriter.open_output_log_for_write(continue_sim /*open in append mode */);
        logwriter.prepare_conservation_file(continue_sim);
        logwriter.prepare_diagnostics_file(continue_sim);
    }


    //*****************************************************************
    //  Set the GPU device.
    cudaError_t error;
    printf(" Using GPU #%d\n", GPU_ID_N);
    cudaSetDevice(GPU_ID_N);

    //
    //  PRINTS
    //  Device Information
    int         ndevices;
    cudaError_t err = cudaGetDeviceCount(&ndevices);

    // Check device query
    if (err != cudaSuccess) {
        printf("Error getting device count.\n");
        printf("%s\n", cudaGetErrorString(err));
        exit(-1);
    }

    int device_major_minor_number = 0;
    for (int i = 0; i < ndevices; ++i) {
        // Get device properties
        printf("\n CUDA Device #%d\n", i);
        cudaDeviceProp devPp;
        cudaGetDeviceProperties(&devPp, i);

        if (i == GPU_ID_N)
            device_major_minor_number = devPp.major * 10 + devPp.minor;

        printf(" Name: %s\n", devPp.name);
        printf(" Compute Capabilities: %d.%d\n", devPp.major, devPp.minor);
        printf("   Total global memory:           %lu\n", devPp.totalGlobalMem);
        printf("   Total shared memory per block: %lu\n", devPp.sharedMemPerBlock);
        printf("   Total registers per block:     %d\n", devPp.regsPerBlock);
        printf("   Warp size:                     %d\n", devPp.warpSize);
        printf("   Maximum memory pitch:          %lu\n", devPp.memPitch);
        printf("   Maximum threads per block:     %d\n", devPp.maxThreadsPerBlock);
        printf("   Clock rate:                    %d\n", devPp.clockRate);
        printf("   Total constant memory:         %lu\n", devPp.totalConstMem);
        printf("   Number of multiprocessors:     %d\n", devPp.multiProcessorCount);
    }

    // do we have a device?
    if (ndevices < 1 || err != cudaSuccess) {
        printf("No device found (compiled SM:%d).\n", DEVICE_SM);
        printf("Aborting.\n");
        exit(-1);
    }

    // can we match the device ID asked?
    if (GPU_ID_N >= ndevices) {
        printf("Asked for device #%d but only found %d devices.\n", GPU_ID_N, ndevices);
        exit(-1);
    }

    // do we have the compute capabilities set at compile time

    if (device_major_minor_number < DEVICE_SM) {
        printf("Found device with id %d does not have sufficent compute capabilities.\n", GPU_ID_N);
        printf("Capabilities: %d (compiled with SM=%d).\n", device_major_minor_number, DEVICE_SM);
        printf("Aborting.\n");
        exit(-1);
    }
    else if (device_major_minor_number > DEVICE_SM) {
        printf("Device has higher compute capability than used at compile time.\n");
        printf("Capabilities: %d (compiled with SM=%d).\n", device_major_minor_number, DEVICE_SM);
    }


    int max_count = 0;
    //
    //  Make the icosahedral grid
    Icogrid Grid(spring_dynamics,     // Spring dynamics option
                 spring_beta,         // Parameter beta for spring dynamics
                 glevel,              // Horizontal resolution level
                 vlevel,              // Number of vertical layers
                 nlat,                // Number of lat rings for sponge layer
                 Planet.A,            // Planet radius
                 Planet.Top_altitude, // Top of the model's domain
                 SpongeLayer,         // Use sponge layer?
                 &max_count);

    //  Define object X.
    ESP X(Grid.point_local,    // First neighbours
          Grid.maps,           // Grid domains
          Grid.lonlat,         // Longitude and latitude of the grid points
          Grid.Altitude,       // Altitudes
          Grid.Altitudeh,      // Altitude at the interfaces between layers
          Grid.nvecoa,         // Normal vectors for diffusion 1
          Grid.nvecti,         // Normal vectors for diffusion 2
          Grid.nvecte,         // Normal vectors for diffusion 3
          Grid.areasT,         // Areas of the main cells
          Grid.areasTr,        // Areas of the triangles
          Grid.div,            // Divergence operator
          Grid.grad,           // Gradient operator
          Grid.func_r,         // Normalised vector
          Grid.nl_region,      // Number of points in one side of a rhombus
          Grid.nr,             // Number of rhombi
          Grid.nv,             // Number of vertical layers
          Grid.nvi,            // Number of interfaces between layer
          glevel,              // Horizontal resolution level
          spring_dynamics,     // Spring dynamics option
          spring_beta,         // Parameter beta for spring dynamics
          nlat,                // Number of latitude rings for zonal
                               // mean wind
          Grid.zonal_mean_tab, // table of zonal means for sponge layer
          Rv_sponge,           // Maximum damping of sponge layer
          ns_sponge,           // lowest level of sponge layer (fraction of model)
          t_shrink,            // time to shrink sponge layer
          Grid.point_num,      // Number of grid points
          conservation,        // compute conservation values
          core_benchmark,      // benchmark test type
          logwriter);          // Log writer

    USE_BENCHMARK();

    INIT_BENCHMARK(X, Grid, output_path_ref);

    BENCH_POINT("0", "Grid", (), ("func_r", "areas", "areasTr", "areasT", "nvec", "nvecoa", "nvecti", "nvecte", "Altitude", "Altitudeh", "lonlat", "div", "grad"))

    // esp output setup
    X.set_output_param(Planet.simulation_ID, output_path);

    printf(" Setting the initial conditions.\n\n");

    double simulation_start_time = 0.0;

    // Initial conditions
    int  output_file_idx = 0;
    int  step_idx        = 0;
    bool load_initial    = X.initial_values(rest,                  // Option to
                                                                // start the
                                                                // atmosphere
                                                                // from rest
                                         initial_conditions,    // initial conditions if not
                                                                // started from
                                                                // rest
                                         continue_sim,          // if we
                                                                // specify
                                                                // initial
                                                                // conditions,
                                                                // continue or
                                                                // start at 0?
                                         timestep,              // Time-step [s]
                                         Planet,                // Planet
                                         SpongeLayer,           // Enable sponge layer
                                         DeepModel,             // Use deep model corrections
                                         TPprof,                // isothermal = 0, guillot = 1
                                         step_idx,              // current step index
                                         simulation_start_time, // output:
                                                                // simulation start time
                                         output_file_idx,       // output file
                                                                // read + 1, 0
                                                                // if nothing read
                                         conservation);

    if (!load_initial) {
        printf("error loading initial conditions from %s.\n", initial_conditions.c_str());
        exit(EXIT_FAILURE);
    }

    // Check presence of output files
    path results(output_path);

    // Get the files in the directory
    std::vector<string> result_files = get_files_in_directory(results.to_string());

    // match them to name pattern, get file numbers and check numbers that are greater than the
    // restart file number
    std::vector<std::pair<string, int>> matching_name_result_files;
    for (const auto& r : result_files) {
        string basename    = "";
        int    file_number = 0;
        if (match_output_file_numbering_scheme(r, basename, file_number)) {
            if (basename == simulation_ID && file_number > output_file_idx)
                matching_name_result_files.emplace_back(r, file_number);
        }
    }
    // sort them by number
    std::sort(matching_name_result_files.begin(),
              matching_name_result_files.end(),
              [](const std::pair<string, int>& left, const std::pair<string, int>& right) {
                  return left.second < right.second;
              });

    if (matching_name_result_files.size() > 0) {
        if (!force_overwrite) {
            printf("output files already exist and would be overwritten \n"
                   "when running simulation. \n"
                   "Files found:\n");
            for (const auto& f : matching_name_result_files)
                printf("\t%s\n", f.first.c_str());

            printf(" Aborting. \n"
                   "use --overwrite to overwrite existing files.\n");
            exit(EXIT_FAILURE);
        }
    }


    long startTime = clock();


    //
    //  Planet conditions
    printf("\n");
    printf(" Planet: %s\n", Planet.simulation_ID);
    printf("   Radius = %f m\n", Planet.A);
    printf("   Omega  = %f s-1\n", Planet.Omega);
    printf("   Gravit = %f m/s2\n", Planet.Gravit);
    printf("   Rd     = %f J/(Kg K)\n", Planet.Rd);
    printf("   Cp     = %f J/(Kg K)\n", Planet.Cp);
    printf("   Tmean  = %f K\n", Planet.Tmean);
    //
    //  Numerical Methods
    printf("\n");
    printf(" THOR\n");
    printf("   ********** \n");
    printf("   Grid - Icosahedral\n");
    printf("   Glevel          = %d.\n", glevel);
    printf("   Spring dynamics = %d.\n", spring_dynamics);
    printf("   Beta            = %f.\n", spring_beta);
    printf("   Resolution      = %f deg.\n", (180 / M_PI) * sqrt(2 * M_PI / 5) / pow(2, glevel));
    printf("   Vertical layers = %d.\n", Grid.nv);
    printf("   ********** \n");
    printf("   Split-Explicit / HE-VI \n");
    printf("   FV = Central finite volume \n");
    printf("   Time integration =  %d s.\n", nsmax * timestep);
    printf("   Large time-step  =  %d s.\n", timestep);
    printf("   Start time       =  %f s.\n", simulation_start_time);

    printf("    \n");

    printf("   Running Core Benchmark test \"%s\" (%d).\n", core_benchmark_str.c_str(), int(core_benchmark));

    printf("    \n");

    printf("\n");
    printf(" Physics module: %s   \n", phy_modules_get_name().c_str());
    printf("   ********** \n");

    if (core_benchmark == NO_BENCHMARK) {

        printf("    \n");
        phy_modules_print_config();
        printf("    \n");
    }

    printf("   ********** \n");
    printf("\n");
    printf("    \n");
    printf(" Simulation\n");

    printf("   Start from rest = %s \n", rest ? "true" : "false");
    if (!rest)
        printf("   Loading initial conditions from = %s \n", initial_conditions.c_str());
    printf("   Output directory = %s \n", output_path.c_str());
    printf("   Start output numbering at %d.\n", output_file_idx);


    // We'll start writnig data to file and running main loop,
    // setup signal handlers to handle gracefully termination and interrupt
    struct sigaction sigterm_action;

    sigterm_action.sa_handler = &sigterm_handler;
    sigemptyset(&sigterm_action.sa_mask);
    sigterm_action.sa_flags = 0;

    if (sigaction(SIGTERM, &sigterm_action, NULL) == -1) {
        printf("Error: cannot handle SIGTERM\n"); // Should not happen
    }

    struct sigaction sigint_action;

    sigint_action.sa_handler = &sigint_handler;
    sigemptyset(&sigint_action.sa_mask);
    sigint_action.sa_flags = 0;

    if (sigaction(SIGINT, &sigint_action, NULL) == -1) {
        printf("Error: cannot handle SIGINT\n"); // Should not happen
    }

    //
    //  Writes initial conditions
    double simulation_time = simulation_start_time;
    if (!continue_sim) {
        X.copy_to_host();
        X.init_timestep(0, simulation_time, timestep);

        if (conservation == true) {
            X.conservation(Planet,
                           DeepModel);

            logwriter.output_conservation(0,
                                          simulation_time,
                                          X.GlobalE_h,
                                          X.GlobalMass_h,
                                          X.GlobalAMx_h,
                                          X.GlobalAMy_h,
                                          X.GlobalAMz_h);
        }

        X.output(0, // file index
                 Planet,
                 conservation,
                 SpongeLayer);
        output_file_idx = 1;
        step_idx        = 1;
    }
    else {
        output_file_idx += 1;
        step_idx += 1;
    }

    // *********************************************************************************************
    //  Starting model Integration.
    printf(" Starting the model integration.\n\n");

    // Start timer
    iteration_timer ittimer(step_idx, nsmax);

    //
    //  Main loop. nstep is the current step of the integration and nsmax the maximum
    //  number of steps in the integration.
    for (int nstep = step_idx; nstep <= nsmax; ++nstep) {

        // compute simulation time
        simulation_time = simulation_start_time + (nstep - step_idx + 1) * timestep;
        // set simulation time and step number for simulation engine and output
        X.init_timestep(nstep,           // Time-step [s]
                        simulation_time, // Simulation time [s]
                        timestep);       // Large time step [s]

        if (!gcm_off) {
            //
            //        Dynamical Core Integration (THOR)
            X.Thor(Planet,     // planet parameters
                   HyDiff,     // Hyperdiffusion option
                   DivDampP,   // Divergence-damping option
                   NonHydro,   // Non-hydrostatic option
                   DeepModel); // Deep model option
        }
        //
        //     Physical Core Integration (ProfX)
        X.ProfX(Planet,
                conv_adj,
                DeepModel,
                n_out,
                SpongeLayer,
                shrink_sponge,
                conservation);

        // compute simulation time
        simulation_time  = simulation_start_time + (nstep - step_idx + 1) * timestep;
        bool file_output = false;

        if (conservation == true) {
            X.conservation(Planet,
                           DeepModel);
            logwriter.output_conservation(nstep,
                                          simulation_time,
                                          X.GlobalE_h,
                                          X.GlobalMass_h,
                                          X.GlobalAMx_h,
                                          X.GlobalAMy_h,
                                          X.GlobalAMz_h);
        }

        //
        //      Prints output every nout steps
        if (nstep % n_out == 0
            || caught_signal != ESIG_NOSIG) {
            X.copy_to_host();
            X.output(output_file_idx,
                     Planet,
                     conservation,
                     SpongeLayer);
            // increment output file index
            output_file_idx++;

            file_output = true;
        }


        // Timing information
        double      mean_delta_per_step = 0.0;
        double      elapsed_time        = 0.0;
        double      time_left           = 0.0;
        std::time_t end_time;

        ittimer.iteration(nstep,
                          mean_delta_per_step,
                          elapsed_time,
                          time_left,
                          end_time);
        // format end time
        std::ostringstream end_time_str;
        char               str_time[256];
        std::strftime(str_time, sizeof(str_time), "%F %T", std::localtime(&end_time));
        end_time_str << str_time;


        printf("\n Time step number = %d/%d || Time = %f days. \n\t Elapsed %s || Left: %s || Completion: %s. %s",
               nstep,
               nsmax,
               simulation_time / 86400.,
               duration_to_str(elapsed_time).c_str(),
               duration_to_str(time_left).c_str(),
               end_time_str.str().c_str(),
               file_output ? "[saved output]" : "");

        // get memory statistics
        size_t total_bytes;
        size_t free_bytes;
        get_cuda_mem_usage(total_bytes, free_bytes);


        logwriter.output_diagnostics(nstep,
                                     simulation_time,
                                     total_bytes,
                                     free_bytes,
                                     elapsed_time,
                                     time_left,
                                     mean_delta_per_step,
                                     end_time);

        if (caught_signal != ESIG_NOSIG) {
            //exit loop and application after save on SIGTERM or SIGINT
            break;
        }
    }
    //
    //  Prints the duration of the integration.
    long finishTime = clock();
    printf("\n\n Integration time = %f seconds\n\n", double((finishTime - startTime)) / double(CLOCKS_PER_SEC));

    //
    //  Checks for errors in the device.
    cudaDeviceSynchronize();
    error = cudaGetLastError();
    printf("CudaMalloc error = %d = %s\n\n", error, cudaGetErrorString(error));
    //
    //  END OF THE ESP
    printf("End of the ESP!\n\n");

    return 0;
}
