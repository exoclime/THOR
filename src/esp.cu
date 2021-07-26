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
// ESP -  Exoclimes Simulation Platform. (version 2.3)
//
//
//
// Method: [1] - Builds the icosahedral grid (Mendonca et al., 2016).
//         Iteration:
//         [2] - Calls the physical core ProfX.
//         [3] - Calls the dynamical core THOR (Mendonca et al., 2016).
//         END
////
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

#include "define.h"
#include "esp.h"
#include "grid.h"
#include "simulation_setup.h"

#include "cmdargs.h"
#include "config_file.h"
#include "directories.h"
#include "iteration_timer.h"

#include "binary_test.h"

#include "phy_modules.h"

#include "debug.h"
#include "diagnostics.h"
#ifdef BENCHMARKING
#    warning "Compiling with benchmarktest enabled"
#endif

#include <csignal>

using std::string;

#include "log_writer.h"

#include "cuda_device_memory.h"

#include "reduction_min.h"

#include "insolation.h"

#include "git-rev.h"

#include "hdf5.h"

enum e_sig { ESIG_NOSIG = 0, ESIG_SIGTERM = 1, ESIG_SIGINT = 2 };


volatile sig_atomic_t caught_signal = ESIG_NOSIG;


void sigterm_handler(int sig) {
    caught_signal = ESIG_SIGTERM;
    log::printf("SIGTERM caught, trying to exit gracefully.\n");
}

void sigint_handler(int sig) {
    caught_signal = ESIG_SIGINT;
    log::printf("SIGINT caught, trying to exit gracefully\n");
}

// To clean up at exit in all conditions (like when calling exit(EXIT_FAILURE)
void exit_handler() {
    // Clean up device memory
    cuda_device_memory_manager::get_instance().deallocate();
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
        log::printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
    }
}


void print_hdf_5_version() {
    unsigned int majnum = 0;
    unsigned int minnum = 0;
    unsigned int relnum = 0;

    H5get_libversion(&majnum, &minnum, &relnum);

    printf("using HDF5 version %d.%d.%d\n", majnum, minnum, relnum);
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

    // for bool, the default value given is the value that will be returned if the option is present.

    string default_config_filename("ifile/earth.thr");

    argparser.add_positional_arg(default_config_filename, "config filename");
    argparser.add_arg("g", "gpu_id", 0, "GPU_ID to run on");
    argparser.add_arg("o", "output_dir", string("results"), "results directory to store output");

    argparser.add_arg("i",
                      "initial",
                      string("initialfilename"),
                      "start from this initial condition instead of rest");
    argparser.add_arg(
        "c", "continue", string("continuefilename"), "continue simulation from this output file");
    argparser.add_arg("N", "numsteps", 48000, "number of steps to run");
    argparser.add_arg("w", "overwrite", true, "Force overwrite of output file if they exist");

    argparser.add_arg("b", "batch", true, "Run as batch");

#ifdef BENCHMARKING
    // binary comparison type
    argparser.add_arg("bw", "binwrite", true, "binary comparison write");
    argparser.add_arg("bc", "bincompare", true, "binary comparison compare");
#endif // BENCHMARKING

    // Parse arguments, exit on fail
    bool parse_ok = argparser.parse(argc, argv);
    if (!parse_ok) {
        log::printf("error parsing arguments\n");

        argparser.print_help();
        exit(-1);
    }


    // get config file path if present in argument
    string config_filename = "";
    argparser.get_positional_arg(config_filename);

#ifdef BENCHMARKING
    // binary comparison reading
    bool bincomp_write     = false;
    bool bincomp_write_arg = false;
    if (argparser.get_arg("binwrite", bincomp_write_arg))
        bincomp_write = bincomp_write_arg;

    bool bincomp_compare     = false;
    bool bincomp_compare_arg = false;
    if (argparser.get_arg("bincompare", bincomp_compare_arg))
        bincomp_compare = bincomp_compare_arg;

    if (bincomp_compare && bincomp_write) {
        log::printf(
            "ARGUMENT ERROR: binwrite and bincompare can't be set to true at the same time\n");
        exit(-1);
    }
#endif // BENCHMARKING

    //*****************************************************************
    log::printf("\n Starting ESP!");
    log::printf("\n\n version 2.3\n");
    log::printf(" Compiled on %s at %s.\n", __DATE__, __TIME__);

    log::printf(" build level: %s\n", BUILD_LEVEL);

    log::printf(" git branch: %s\n", GIT_BRANCH_RAW);
    log::printf(" git hash: %s\n", GIT_HASH_RAW);


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

    // initialises to default simulation, default Earth parameters
    SimulationSetup sim;

    string simulation_ID = "Earth";

    config_reader.append_config_var("simulation_ID", simulation_ID, string(sim.simulation_ID));
    config_reader.append_config_var("radius", sim.A, sim.A);
    config_reader.append_config_var("rotation_rate", sim.Omega, sim.Omega);
    config_reader.append_config_var("gravitation", sim.Gravit, sim.Gravit);
    config_reader.append_config_var("Rd", sim.Rd, sim.Rd);
    config_reader.append_config_var("Cp", sim.Cp, sim.Cp);
    config_reader.append_config_var("Tmean", sim.Tmean, sim.Tmean);
    config_reader.append_config_var("P_ref", sim.P_Ref, sim.P_Ref);
    config_reader.append_config_var("Top_altitude", sim.Top_altitude, sim.Top_altitude);
    config_reader.append_config_var("Diffc", sim.Diffc, sim.Diffc);

    // ...and separate div damping
    config_reader.append_config_var("DivDampc", sim.DivDampc, sim.DivDampc);
    // vertical hyper diffusion
    config_reader.append_config_var("Diffc_v", sim.Diffc_v, sim.Diffc_v);


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

    config_reader.append_config_var("HyDiff", sim.HyDiff, HyDiff_default);
    config_reader.append_config_var("HyDiffOrder", sim.HyDiffOrder, HyDiffOrder_default);
    config_reader.append_config_var("DivDampP", sim.DivDampP, DivDampP_default);
    config_reader.append_config_var("VertHyDiff", sim.VertHyDiff, VertHyDiff_default);
    config_reader.append_config_var(
        "VertHyDiffOrder", sim.VertHyDiffOrder, VertHyDiffOrder_default);


    // Model options
    config_reader.append_config_var("NonHydro", sim.NonHydro, NonHydro_default);
    config_reader.append_config_var("DeepModel", sim.DeepModel, DeepModel_default);
    config_reader.append_config_var("output_mean", sim.output_mean, output_mean_default);
    config_reader.append_config_var(
        "output_diffusion", sim.output_diffusion, output_diffusion_default);
    config_reader.append_config_var(
        "out_interm_momentum", sim.out_interm_momentum, out_interm_momentum_default);

    // top sponge layer options
    int    nlat_bins;
    double Ruv_sponge;
    double Rw_sponge;
    double RT_sponge;
    double ns_ray_sponge;
    bool   shrink_sponge;
    int    t_shrink; // number of time steps after which shrink begins
    bool   damp_uv_to_mean;
    bool   damp_w_to_mean;
    double Dv_sponge;
    double ns_diff_sponge;
    int    order_diff_sponge;
    string raysp_calc_mode_str("imp");

    config_reader.append_config_var("RayleighSponge", sim.RayleighSponge, RayleighSponge_default);
    config_reader.append_config_var(
        "RayleighSpongeT", sim.RayleighSpongeT, RayleighSpongeT_default);
    config_reader.append_config_var("DiffSponge", sim.DiffSponge, DiffSponge_default);
    config_reader.append_config_var("shrink_sponge", shrink_sponge, shrink_sponge_default);
    config_reader.append_config_var("t_shrink", t_shrink, t_shrink_default);

    config_reader.append_config_var(
        "raysp_calc_mode", raysp_calc_mode_str, string(raysp_calc_mode_default));
    config_reader.append_config_var("nlat_bins", nlat_bins, nlat_bins_default);
    config_reader.append_config_var("Ruv_sponge", Ruv_sponge, Ruv_sponge_default);
    config_reader.append_config_var("Rw_sponge", Rw_sponge, Rw_sponge_default);
    config_reader.append_config_var("RT_sponge", RT_sponge, RT_sponge_default);
    config_reader.append_config_var("ns_ray_sponge", ns_ray_sponge, ns_ray_sponge_default);
    config_reader.append_config_var("damp_uv_to_mean", damp_uv_to_mean, damp_uv_to_mean_default);
    config_reader.append_config_var("damp_w_to_mean", damp_w_to_mean, damp_w_to_mean_default);
    config_reader.append_config_var("Dv_sponge", Dv_sponge, Dv_sponge_default);
    config_reader.append_config_var("ns_diff_sponge", ns_diff_sponge, ns_diff_sponge_default);
    config_reader.append_config_var(
        "order_diff_sponge", order_diff_sponge, order_diff_sponge_default);

    // Low pressure test
    bool   exit_on_low_pressure_warning = false;
    double pressure_check_limit         = 1e-3;
    config_reader.append_config_var(
        "exit_on_low_pressure_warning", exit_on_low_pressure_warning, exit_on_low_pressure_warning);
    config_reader.append_config_var(
        "pressure_check_limit", pressure_check_limit, pressure_check_limit);
    // Initial conditions
    // rest supersedes initial condition entry,
    // but if initial condition set from  command line, it overrides
    // rest variable
    config_reader.append_config_var("rest", sim.rest, rest_default);

    string initial_conditions = "initialfilename.h5";
    config_reader.append_config_var(
        "initial", initial_conditions, string(initial_conditions_default));

    // Benchmark test
    string core_benchmark_str("HeldSuarez");
    config_reader.append_config_var(
        "core_benchmark", core_benchmark_str, string(core_benchmark_default)); //

    config_reader.append_config_var("conv_adj", sim.conv_adj, conv_adj_default);
    config_reader.append_config_var("conv_adj_iter", sim.conv_adj_iter, conv_adj_iter_default);

    int GPU_ID_N = 0;
    config_reader.append_config_var("GPU_ID_N", GPU_ID_N, GPU_ID_N_default);

    // int n_out = 1000;
    config_reader.append_config_var("n_out", sim.n_out, n_out_default);

    string output_path = "results";
    config_reader.append_config_var("results_path", output_path, string(output_path_default));

    config_reader.append_config_var("gcm_off", sim.gcm_off, gcm_off_default);
    config_reader.append_config_var("globdiag", sim.globdiag, globdiag_default);

    config_reader.append_config_var("single_column", sim.single_column, single_column_default);

    bool custom_global_n_out;
    config_reader.append_config_var(
        "custom_global_n_out", custom_global_n_out, custom_global_n_out_default);
    int global_n_out = sim.n_out;
    config_reader.append_config_var("global_n_out", global_n_out, global_n_out_default);

    bool custom_log_n_out;
    config_reader.append_config_var("custom_log_n_out", custom_log_n_out, custom_log_n_out_default);
    int log_n_out = sim.n_out;
    config_reader.append_config_var("log_n_out", log_n_out, log_n_out_default);

    // Init PT profile
    string init_PT_profile_str("isothermal");
    config_reader.append_config_var(
        "init_PT_profile", init_PT_profile_str, string(init_PT_profile_default)); //
    // additional settings for guillot profile (also borrowed by double gray RT)
    double Tint = 100, kappa_lw = 0.002, kappa_sw = 0.001, f_lw = 0.5, bv_freq = 0.01;
    config_reader.append_config_var("Tint", Tint, Tint_default);
    config_reader.append_config_var("kappa_lw", kappa_lw, kappa_lw_default);
    config_reader.append_config_var("kappa_sw", kappa_sw, kappa_sw_default);
    config_reader.append_config_var("f_lw", f_lw, f_lw_default);
    config_reader.append_config_var("bv_freq", bv_freq, bv_freq_default);

    // ultrahot thermodynamics
    string uh_thermo_str("none");
    config_reader.append_config_var("ultrahot_thermo", uh_thermo_str, string(uh_thermo_default)); //

    // ultrahot heating from H-H2 dissociation/recombination
    string uh_heating_str("none");
    config_reader.append_config_var(
        "ultrahot_heating", uh_heating_str, string(uh_heating_default)); //

    string thermo_equation_str("entropy");
    config_reader.append_config_var(
        "thermo_equation", thermo_equation_str, string(thermo_equation_default));

    // vertical grid refinement in lower atmos (for turbulent BL)
    bool vert_refined = false;
    //int  n_bl_layers  = 9;
    config_reader.append_config_var("vert_refined", vert_refined, vert_refined_default);
    //config_reader.append_config_var("n_bl_layers", n_bl_layers, n_bl_layers_default);
    double transition_altitude = 1000; //height of transition from exponential to linear spacing
    double lowest_layer_thickness =
        2; //height of first interface above surface (thickness of lowest layer)
    config_reader.append_config_var(
        "transition_altitude", transition_altitude, transition_altitude_default);
    config_reader.append_config_var(
        "lowest_layer_thickness", lowest_layer_thickness, lowest_layer_thickness_default);


    bool surface_config = false; // use solid/liquid surface at altitude 0
    config_reader.append_config_var("surface", surface_config, surface_config);
    // properties for a solid/liquid surface
    double Csurf_config = 1e7; // heat capacity of surface (J K^-1 m^-2)
    config_reader.append_config_var("Csurf", Csurf_config, Csurf_config);


    //*****************************************************************
    // set configs for modules
    phy_modules_generate_config(config_reader);

    // Create insolation object
    Insolation insolation;
    // set config for insolation
    // insolation will be enabled on demand by modules needing it.
    insolation.configure(config_reader);

    //*****************************************************************
    // Read config file
    log::printf("\n");

    if (config_reader.parse_file(config_filename))
        log::printf(" Config file %s read\n", config_filename.c_str());
    else {
        log::printf(" Config file %s reading failed, aborting\n", config_filename.c_str());
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
        sim.rest                  = false;
        initial_conditions        = inital_conditions_arg;
        initial_condition_arg_set = true;

        if (!path_exists(initial_conditions)) {
            log::printf("Initial conditions file \"%s\" not found\n", initial_conditions.c_str());
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
        sim.rest     = false;
        continue_sim = true;

        if (run_as_batch) {
            log::printf(
                "--continue and --batch options set, options are exclusive, must set only one\n");

            exit(-1);
        }


        if (initial_condition_arg_set) {
            log::printf(
                "--continue and --initial options set, options are exclusive, must set only one\n");

            exit(-1);
        }

        if (!path_exists(continue_filename)) {
            log::printf("Continuation start condition file \"%s\" not found\n",
                        continue_filename.c_str());
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
    config_OK &= check_greater("gravitation", sim.Gravit, 0.0);
    config_OK &= check_greater("Rd", sim.Rd, 0.0);
    config_OK &= check_greater("T_mean", sim.Tmean, 0.0);
    config_OK &= check_greater("P_Ref", sim.P_Ref, 0.0);
    config_OK &= check_greater("Top_altitude", sim.Top_altitude, 0.0);

    config_OK &= check_range("glevel", glevel, 3, 8);
    config_OK &= check_greater("vlevel", vlevel, 0);

    config_OK &= check_greater("GPU_ID_N", GPU_ID_N, -1);
    config_OK &= check_greater("n_out", sim.n_out, 0);

    if (simulation_ID.length() < 160) {
        sprintf(sim.simulation_ID, "%s", simulation_ID.c_str());
    }
    else {
        log::printf("Bad value for config variable simulation_ID: [%s]\n", simulation_ID.c_str());

        config_OK = false;
    }

    if (sim.HyDiffOrder <= 2 || sim.HyDiffOrder % 2) {
        log::printf("Bad value for HyDiffOrder! Must be multiple of 2 greater than/equal to 4.");
        config_OK = false;
    }

    if (sim.VertHyDiffOrder <= 2 || sim.VertHyDiffOrder % 2) {
        log::printf(
            "Bad value for VertHyDiffOrder! Must be multiple of 2 greater than/equal to 4.");
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
        log::printf("core_benchmark 'JetSteady' not fully implemented yet\n");
        config_OK &= false;
    }
    else if (core_benchmark_str == "AcousticTest") {
        core_benchmark = ACOUSTIC_TEST;
        config_OK &= true;
    }
    else if (core_benchmark_str == "GWaveTest") {
        core_benchmark = GWAVE_TEST;
        config_OK &= true;
    }
    else {
        log::printf("core_benchmark config item not recognised: [%s]\n",
                    core_benchmark_str.c_str());
        config_OK &= false;
    }

    // check init pt profile string
    init_PT_profile_types init_PT_profile = ISOTHERMAL;

    if (init_PT_profile_str == "isothermal") {
        init_PT_profile = ISOTHERMAL;
        config_OK &= true;
    }
    
    else if (init_PT_profile_str == "parmentier") {
        init_PT_profile = PARMENTIER;
        config_OK &= true;
    }
    
    else if (init_PT_profile_str == "guillot") {
        init_PT_profile = GUILLOT;
        config_OK &= true;
    }
    else if (init_PT_profile_str == "constbv") {
        init_PT_profile = CONSTBV;
        config_OK &= true;
    }
    else {
        log::printf("init_PT_profile config item not recognised: [%s]\n",
                    init_PT_profile_str.c_str());
        config_OK &= false;
    }

    thermo_equation_types thermo_equation = ENTROPY;

    if (thermo_equation_str == "entropy") {
        thermo_equation = ENTROPY;
        config_OK &= true;
    }
    else if (thermo_equation_str == "energy") {
        thermo_equation = ENERGY;
        config_OK &= true;
    }
    else {
        log::printf("thermo_equation config item not recognised: [%s]\n",
                    thermo_equation_str.c_str());
        config_OK &= false;
    }

    // check ultrahot thermo string
    uh_thermo_types ultrahot_thermo = NO_UH_THERMO;

    if (uh_thermo_str == "none") {
        ultrahot_thermo = NO_UH_THERMO;
        config_OK &= true;
    }
    else if (uh_thermo_str == "vary_R_CP") {
        ultrahot_thermo = VARY_R_CP;
        log::printf("ultrahot_thermo option 'vary_R_CP' not ready yet \n");
        config_OK &= false;
    }
    else if (uh_thermo_str == "full") {
        ultrahot_thermo = FULL;
        log::printf("ultrahot_thermo option 'full' not ready yet \n");
        config_OK &= false;
    }
    else {
        log::printf("ultrahot_thermo config item not recognised: [%s]\n", uh_thermo_str.c_str());
        config_OK &= false;
    }

    // check ultrahot heating string
    uh_heating_types ultrahot_heating = NO_UH_HEATING;

    if (uh_heating_str == "none") {
        ultrahot_heating = NO_UH_HEATING;
        config_OK &= true;
    }
    else if (uh_heating_str == "quasi_eql") {
        ultrahot_heating = QUASI_EQL;
        log::printf("ultrahot_heating option 'quasi_eql' not ready yet \n");
        config_OK &= false;
    }
    else if (uh_heating_str == "relax_chem") {
        ultrahot_heating = RELAX_CHEM;
        log::printf("ultrahot_heating option 'relax_chem' not ready yet \n");
        config_OK &= false;
    }
    else {
        log::printf("ultrahot_heating config item not recognised: [%s]\n", uh_heating_str.c_str());
        config_OK &= false;
    }

    bool                  initialize_zonal_mean = false;
    raysp_calc_mode_types raysp_calc_mode       = IMP;
    if (sim.RayleighSponge) {

        if (raysp_calc_mode_str == "imp") {
            raysp_calc_mode = IMP;
            config_OK &= true;
        }
        else if (raysp_calc_mode_str == "exp1") {
            raysp_calc_mode = EXP1;
            config_OK &= true;
        }
        else if (raysp_calc_mode_str == "exp3") {
            raysp_calc_mode = EXP3;
            config_OK &= true;
        }
        else {
            log::printf("raysp_calc_mode config item not recognised: [%s]\n",
                        raysp_calc_mode_str.c_str());
            config_OK &= false;
        }

        if (damp_uv_to_mean || damp_w_to_mean) {
            initialize_zonal_mean = true;
        }
    }

    if (sim.DiffSponge) {
        if (order_diff_sponge == 2 || order_diff_sponge == 4) {
            config_OK &= true;
        }
        else {
            log::printf("order_diff_sponge config option can only be 2 or 4\n");
            config_OK &= false;
        }
    }

    if ((sim.single_column == true) && (sim.gcm_off == false)) {
        log::printf("gcm_off must be true when using single_column\n");
        config_OK &= false;
    }

    if (!config_OK) {
        log::printf("Error in configuration file\n");
        exit(-1);
    }


#ifdef BENCHMARKING
    string output_path_ref = path(output_path).to_string();
    if (bincomp_compare) {
        log::printf("\nWARNING: Running with binary comparison ON\n\n");
        output_path = (path(output_path) / string("compare")).to_string();
    }
    if (bincomp_write) {
        log::printf("\nWARNING: Running with binary comparison data write ON\n\n");
        output_path = (path(output_path) / string("write")).to_string();
    }
#endif // BENCHMARKING


    // check output config directory
    if (!create_output_dir(output_path)) {
        log::printf("Error creating output result directory: %s\n", output_path.c_str());
        exit(-1);
    }

    //****************************************************************
    // Start logging to file
    printf("Opening log file.\n");

    string log_file =
        (path(output_path) / (string("esp_log_") + sim.simulation_ID + string(".log"))).to_string();

    log::init_logger(log_file, continue_sim || run_as_batch);

    log::printf_logonly("*********************************************************************\n");
    std::time_t prog_start_time = std::time(nullptr);
    log::printf("Starting logging to %s.\n\n", log_file.c_str());
    log::printf("Time: %s\n\n", std::ctime(&prog_start_time));


    log::printf_logonly("\n\n version 2.3\n");
    log::printf_logonly(" Compiled on %s at %s.\n", __DATE__, __TIME__);

    log::printf_logonly(" build level: %s\n\n", BUILD_LEVEL);
    log::printf_logonly(" git branch: %s\n", GIT_BRANCH_RAW);
    log::printf_logonly(" git hash: %s\n", GIT_HASH_RAW);


    // Data logger
    log_writer logwriter(sim.simulation_ID, output_path);

    string planet_filename;
    int    output_file_idx = 0;

    // Batch mode handling
    if (run_as_batch) {
        log::printf("Starting in batch mode.\n");


        // Get last written file from
        int    last_file_number      = 0;
        int    last_iteration_number = 0;
        string last_file             = "";
        bool   has_last_file         = false;

        try {
            has_last_file =
                logwriter.check_output_log(last_file_number, last_iteration_number, last_file);
        } catch (const std::exception& e) {
            log::printf(
                "[%s:%d] error while checking output log: %s.\n", __FILE__, __LINE__, e.what());
            exit(-1);
        }

        if (has_last_file) {
            path o(output_path);
            o /= last_file;

            if (path_exists(o.to_string())) {
                log::printf("continuing batch run from file %s\n", last_file.c_str());

                // we want to continue the simulation
                sim.rest     = false;
                continue_sim = true;

                // overwrite the results files we'd encounter
                force_overwrite = true;

                // reload the last file we found as initial conditions
                initial_conditions = o.to_string();

                bool find_continue_OK = false;
                find_continue_OK      = find_continue_file(
                    initial_conditions, planet_filename, continue_sim, output_file_idx);

                if (!find_continue_OK) {
                    cuda_device_memory_manager::get_instance().deallocate();

                    exit(EXIT_FAILURE);
                }

                logwriter.open_output_log_for_write(true /*open in append mode */);
                logwriter.prepare_globdiag_file(true);
                logwriter.prepare_diagnostics_file(true);
            }
            else {
                log::printf("Did not find last saved file that should exist.\n");
                exit(-1);
            }
        }
        else {
            bool find_continue_OK = false;
            if (sim.rest == false) {
                find_continue_OK = find_continue_file(
                    initial_conditions, planet_filename, continue_sim, output_file_idx);
            }
            else {
                find_continue_OK = true;
            }

            if (!find_continue_OK) {
                cuda_device_memory_manager::get_instance().deallocate();

                exit(EXIT_FAILURE);
            }

            log::printf("No batch file found, initialise simulation.\n");
            // we don't have an output file, start from scratch, reinitialising outputs
            logwriter.open_output_log_for_write(false /*open in non append mode */);
            logwriter.prepare_globdiag_file(false);
            logwriter.prepare_diagnostics_file(false);
        }
    }
    else {
        bool find_continue_OK = false;
        if (sim.rest == false) {
            find_continue_OK = find_continue_file(
                initial_conditions, planet_filename, continue_sim, output_file_idx);
        }
        else {
            find_continue_OK = true;
        }

        bool overwrite_OK = false;
        overwrite_OK =
            overwrite_check(output_path, simulation_ID, output_file_idx, force_overwrite);
        if (!overwrite_OK || !find_continue_OK) {
            cuda_device_memory_manager::get_instance().deallocate();

            exit(EXIT_FAILURE);
        }

        log::printf("Opening result output file.\n");
        logwriter.open_output_log_for_write(continue_sim /*open in append mode */);
        logwriter.prepare_globdiag_file(continue_sim);
        logwriter.prepare_diagnostics_file(continue_sim);
    }

    //*****************************************************************
    // print HDF5 version
    print_hdf_5_version();
    //*****************************************************************
    //  Set the GPU device.
    cudaError_t error;
    log::printf(" Using GPU #%d\n", GPU_ID_N);
    cudaSetDevice(GPU_ID_N);

    //
    //  PRINTS
    //  Device Information
    int         ndevices;
    cudaError_t err = cudaGetDeviceCount(&ndevices);

    // Check device query
    if (err != cudaSuccess) {
        log::printf("Error getting device count.\n");
        log::printf("%s\n", cudaGetErrorString(err));
        exit(-1);
    }

    int device_major_minor_number = 0;
    for (int i = 0; i < ndevices; ++i) {
        // Get device properties
        log::printf("\n CUDA Device #%d\n", i);
        cudaDeviceProp devPp;
        cudaGetDeviceProperties(&devPp, i);

        if (i == GPU_ID_N)
            device_major_minor_number = devPp.major * 10 + devPp.minor;

        log::printf(" Name: %s\n", devPp.name);
        log::printf(" Compute Capabilities: %d.%d\n", devPp.major, devPp.minor);
        log::printf("   Total global memory:           %lu\n", devPp.totalGlobalMem);
        log::printf("   Total shared memory per block: %lu\n", devPp.sharedMemPerBlock);
        log::printf("   Total registers per block:     %d\n", devPp.regsPerBlock);
        log::printf("   Warp size:                     %d\n", devPp.warpSize);
        log::printf("   Maximum memory pitch:          %lu\n", devPp.memPitch);
        log::printf("   Maximum threads per block:     %d\n", devPp.maxThreadsPerBlock);
        log::printf("   Clock rate:                    %d\n", devPp.clockRate);
        log::printf("   Total constant memory:         %lu\n", devPp.totalConstMem);
        log::printf("   Number of multiprocessors:     %d\n", devPp.multiProcessorCount);
    }

    // do we have a device?
    if (ndevices < 1 || err != cudaSuccess) {
        log::printf("No device found (compiled SM:%d).\n", DEVICE_SM);
        log::printf("Aborting.\n");
        exit(-1);
    }

    // can we match the device ID asked?
    if (GPU_ID_N >= ndevices) {
        log::printf("Asked for device #%d but only found %d devices.\n", GPU_ID_N, ndevices);
        exit(-1);
    }

    // do we have the compute capabilities set at compile time

    if (device_major_minor_number < DEVICE_SM) {
        log::printf("Found device with id %d does not have sufficent compute capabilities.\n",
                    GPU_ID_N);
        log::printf(
            "Capabilities: %d (compiled with SM=%d).\n", device_major_minor_number, DEVICE_SM);
        log::printf("Aborting.\n");
        exit(-1);
    }
    else if (device_major_minor_number > DEVICE_SM) {
        log::printf("Device has higher compute capability than used at compile time.\n");
        log::printf(
            "Capabilities: %d (compiled with SM=%d).\n", device_major_minor_number, DEVICE_SM);
    }

    int max_count = 0;
    //
    //  Make the icosahedral grid

    Icogrid Grid(spring_dynamics,       // Spring dynamics option
                 spring_beta,           // Parameter beta for spring dynamics
                 glevel,                // Horizontal resolution level
                 vlevel,                // Number of vertical layers
                 nlat_bins,             // Number of lat rings for sponge layer
                 sim.A,                 // Planet radius
                 sim.Top_altitude,      // Top of the model's domain
                 initialize_zonal_mean, // Use zonal mean in rayleigh sponge layer?
                 &max_count,
                 vert_refined,
                 lowest_layer_thickness,
                 transition_altitude);

    //  Define object X.
    int point_num_temp = Grid.point_num;
    if (sim.single_column)
        point_num_temp = 1;

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
          Grid.areas,          // Areas of the sub-triangles
          Grid.div,            // Divergence operator
          Grid.grad,           // Gradient operator
          Grid.curlz,          // Curl operator (vertical component)
          Grid.func_r,         // Normalised vector
          Grid.nl_region,      // Number of points in one side of a rhombus
          Grid.nr,             // Number of rhombi
          Grid.nv,             // Number of vertical layers
          Grid.nvi,            // Number of interfaces between layer
          glevel,              // Horizontal resolution level
          spring_dynamics,     // Spring dynamics option
          spring_beta,         // Parameter beta for spring dynamics
          nlat_bins,           // Number of latitude rings for zonal
                               // mean wind
          Grid.zonal_mean_tab, // table of zonal means for sponge layer
          Ruv_sponge,          // Maximum damping of sponge layer
          Rw_sponge,           // Maximum damping of sponge layer
          RT_sponge,
          ns_ray_sponge, // lowest level of rayleigh sponge layer (fraction of model)
          damp_uv_to_mean,
          damp_w_to_mean,
          raysp_calc_mode,
          Dv_sponge,
          ns_diff_sponge,
          order_diff_sponge,
          t_shrink, // time to shrink sponge layer
          shrink_sponge,
          point_num_temp, // Number of grid points
          sim.globdiag,   // compute globdiag values
          core_benchmark, // benchmark test type
          logwriter,      // Log writer
          max_count,
          sim.output_mean,
          sim.out_interm_momentum,
          sim.output_diffusion,
          init_PT_profile,
          Tint,
          kappa_lw,
          kappa_sw,
          f_lw,
          bv_freq,
          ultrahot_thermo,
          ultrahot_heating,
          thermo_equation,
          surface_config,
          Csurf_config,
          insolation);

    USE_BENCHMARK();

    INIT_BENCHMARK(X, Grid, output_path_ref, bincomp_write, bincomp_compare);

    BENCH_POINT("0",
                "Grid",
                (),
                ("func_r",
                 "areas",
                 "areasTr",
                 "areasT",
                 "nvec",
                 "nvecoa",
                 "nvecti",
                 "nvecte",
                 "Altitude",
                 "Altitudeh",
                 "lonlat",
                 "div",
                 "grad"))

    // esp output setup
    X.set_output_param(sim.simulation_ID, output_path);

    log::printf(" Setting the initial conditions.\n\n");

    double simulation_start_time = 0.0;

    // Initial conditions
    int  step_idx     = 0;
    bool load_initial = X.initial_values(initial_conditions, // initial conditions if not
                                                             // started from
                                                             // rest
                                         planet_filename,
                                         // continue_sim,          // if we
                                         // specify
                                         // initial
                                         // conditions,
                                         // continue or
                                         // start at 0?
                                         timestep,               // Time-step [s]
                                         sim,                    // simulation parameters
                                         step_idx,               // current step index
                                         simulation_start_time); // output:
                                                                 // simulation start time
    // output_file_idx);      // output file
    // read + 1, 0
    // if nothing read


    if (!load_initial) {
        log::printf("error loading initial conditions from %s.\n", initial_conditions.c_str());

        exit(EXIT_FAILURE);
    }

    long startTime = clock();


    //
    //  Planet conditions
    log::printf("\n");
    log::printf(" Planet: %s\n", sim.simulation_ID);
    log::printf("   Radius = %f m\n", sim.A);
    log::printf("   Omega  = %f s-1\n", sim.Omega);
    log::printf("   Gravit = %f m/s2\n", sim.Gravit);
    log::printf("   Rd     = %f J/(Kg K)\n", sim.Rd);
    log::printf("   Cp     = %f J/(Kg K)\n", sim.Cp);
    log::printf("   Tmean  = %f K\n", sim.Tmean);
    //
    //  Numerical Methods
    log::printf("\n");
    log::printf(" THOR\n");
    log::printf("   ********** \n");
    log::printf("   Grid - Icosahedral\n");
    log::printf("   Glevel          = %d.\n", glevel);
    log::printf("   Spring dynamics = %d.\n", spring_dynamics);
    log::printf("   Beta            = %f.\n", spring_beta);
    log::printf("   Resolution      = %f deg.\n",
                (180 / M_PI) * sqrt(2 * M_PI / 5) / pow(2, glevel));
    log::printf("   Vertical layers = %d.\n", Grid.nv);
    log::printf("   ********** \n");
    log::printf("   Split-Explicit / HE-VI \n");
    log::printf("   FV = Central finite volume \n");
    log::printf("   Time integration =  %d s.\n", nsmax * timestep);
    log::printf("   Large time-step  =  %d s.\n", timestep);
    log::printf("   Start time       =  %f s.\n", simulation_start_time);
    // surface parameters
    log::printf("   Surface          = %s.\n", surface_config ? "true" : "false");
    log::printf("   Surface Heat Capacity       = %f J/K/m^2.\n", Csurf_config);

    log::printf("    \n");

    log::printf("   Running Core Benchmark test \"%s\" (%d).\n",
                core_benchmark_str.c_str(),
                int(core_benchmark));

    log::printf("    \n");

    log::printf("\n");
    log::printf(" Physics module: %s   \n", phy_modules_get_name().c_str());
    log::printf("   ********** \n");

    if (core_benchmark == NO_BENCHMARK) {

        log::printf("    \n");
        phy_modules_print_config();
        log::printf("    \n");
        insolation.print_config();
    }

    log::printf("   ********** \n");
    log::printf("\n");
    log::printf("    \n");
    log::printf(" Simulation\n");

    log::printf("   Start from rest = %s \n", sim.rest ? "true" : "false");
    if (!sim.rest)
        log::printf("   Loading initial conditions from = %s \n", initial_conditions.c_str());
    log::printf("   Output directory = %s \n", output_path.c_str());
    log::printf("   Start output numbering at %d.\n", output_file_idx);

    // make a copy of config file (is there a better way to do this?)
    std::ifstream source(config_filename);
    char          dest_name[512];
    sprintf(dest_name, "%s/config_copy.%d", output_path.c_str(), output_file_idx);
    std::ofstream destin(dest_name);

    std::istreambuf_iterator<char> begin_source(source);
    std::istreambuf_iterator<char> end_source;
    std::ostreambuf_iterator<char> begin_dest(destin);
    std::copy(begin_source, end_source, begin_dest);

    source.close();
    destin.close();


    // We'll start writnig data to file and running main loop,
    // setup signal handlers to handle gracefully termination and interrupt
    struct sigaction sigterm_action;

    sigterm_action.sa_handler = &sigterm_handler;
    sigemptyset(&sigterm_action.sa_mask);
    sigterm_action.sa_flags = 0;

    if (sigaction(SIGTERM, &sigterm_action, NULL) == -1) {
        log::printf("Error: cannot handle SIGTERM\n"); // Should not happen
    }

    struct sigaction sigint_action;

    sigint_action.sa_handler = &sigint_handler;
    sigemptyset(&sigint_action.sa_mask);
    sigint_action.sa_flags = 0;

    if (sigaction(SIGINT, &sigint_action, NULL) == -1) {
        log::printf("Error: cannot handle SIGINT\n"); // Should not happen
    }

    // Register clean up function for exit
    atexit(exit_handler);

    //
    //  Writes initial conditions
    double simulation_time = simulation_start_time;
    if (!continue_sim) {
        X.copy_to_host();
        X.init_timestep(0, simulation_time, timestep);

        if (sim.globdiag == true) {
            X.globdiag(sim);


            logwriter.output_globdiag(0,
                                      simulation_time,
                                      X.GlobalE_h,
                                      X.GlobalMass_h,
                                      X.GlobalAMx_h,
                                      X.GlobalAMy_h,
                                      X.GlobalAMz_h,
                                      X.GlobalEnt_h);

            X.copy_globdiag_to_host();
        }

        if (sim.output_mean == true) {
            X.copy_mean_to_host();
        }

        if (sim.out_interm_momentum == true) {
            X.copy_interm_mom_to_host();
        }

        if (sim.output_diffusion == true) {
            X.copy_diff_to_host();
        }

        X.output(0, // file index
                 sim);

        output_file_idx = 1;
        step_idx        = 1;
    }
    else {
        output_file_idx += 1;
        step_idx += 1;
    }

    // Create diagnostics tool
    kernel_diagnostics diag(X);

    // *********************************************************************************************
    //  Starting model Integration.
    log::printf(" Starting the model integration.\n\n");

    // Start timer
    iteration_timer ittimer(step_idx, nsmax);

    //
    //  Main loop. nstep is the current step of the integration and nsmax the maximum
    //  number of steps in the integration.
    int n_since_output = 1; //how many steps have passed since last output
    for (int nstep = step_idx; nstep <= nsmax; ++nstep) {
        // compute simulation time
        simulation_time = simulation_start_time + (nstep - step_idx + 1) * timestep;

        // set simulation time and step number for simulation engine and output
        X.init_timestep(nstep,           // Time-step [s]
                        simulation_time, // Simulation time [s]
                        timestep);       // Large time step [s]

        //
        //     Physical Core Integration (ProfX)
        X.ProfX(sim, sim.n_out, diag);

        if (!sim.gcm_off) {
            //
            //        Dynamical Core Integration (THOR)
            X.Thor(sim, diag); // simulationt parameters
        }

        //
        //     Physical Core Integration (ProfX)
        // X.ProfX(sim, n_out, shrink_sponge);

        // compute simulation time
        simulation_time  = simulation_start_time + (nstep - step_idx + 1) * timestep;
        bool file_output = false;

        if (sim.globdiag == true) {
            if ((custom_global_n_out && nstep % global_n_out == 0)
                || (custom_global_n_out == false && nstep % sim.n_out == 0)) {
                X.globdiag(sim);
                logwriter.output_globdiag(nstep,
                                          simulation_time,
                                          X.GlobalE_h,
                                          X.GlobalMass_h,
                                          X.GlobalAMx_h,
                                          X.GlobalAMy_h,
                                          X.GlobalAMz_h,
                                          X.GlobalEnt_h);
            }
        }

        if (sim.output_mean == true)
            X.update_mean_outputs(n_since_output);

        //
        //      Prints output every nout steps
        if (nstep % sim.n_out == 0 || caught_signal != ESIG_NOSIG) {
            X.copy_to_host();

            if (sim.globdiag == true)
                X.copy_globdiag_to_host();

            if (sim.output_mean == true) {
                X.copy_mean_to_host();
            }

            if (sim.out_interm_momentum == true) {
                X.copy_interm_mom_to_host();
            }

            if (sim.output_diffusion == true) {
                X.copy_diff_to_host();
            }

            X.output(output_file_idx, sim);

            // increment output file index
            output_file_idx++;

            file_output    = true;
            n_since_output = 0;
        }

        // Timing information
        double      mean_delta_per_step = 0.0;
        double      this_step_delta     = 0.0;
        double      elapsed_time        = 0.0;
        double      time_left           = 0.0;
        std::time_t end_time;

        ittimer.iteration(
            nstep, mean_delta_per_step, this_step_delta, elapsed_time, time_left, end_time);
        // format end time
        std::ostringstream end_time_str;
        char               str_time[256];
        std::strftime(str_time, sizeof(str_time), "%F %T", std::localtime(&end_time));
        end_time_str << str_time;

        if ((custom_log_n_out && nstep % log_n_out == 0)
            || (custom_log_n_out == false && nstep % sim.n_out == 0)) {
            log::printf(
                "\n Time step number = %d/%d || Time = %f days. \n\t Elapsed %s || Left: %s || "
                "Completion: %s. %s",
                nstep,
                nsmax,
                simulation_time / 86400.,
                duration_to_str(elapsed_time).c_str(),
                duration_to_str(time_left).c_str(),
                end_time_str.str().c_str(),
                file_output ? "[saved output]" : "");
#ifdef STEP_TIMING_INFO
            log::printf("\n\t Step duration: %f s || Mean step duration %f s",
                        this_step_delta,
                        mean_delta_per_step);
#endif // STEP_TIMING_INFO

            // get memory statistics
            size_t total_bytes;
            size_t free_bytes;
            get_cuda_mem_usage(total_bytes, free_bytes);

            // Run pressure check
            double pressure_min = gpu_min_on_device<1024>(X.pressure_d, X.point_num * X.nv);
            log::printf("\n min pressure : %g Pa\n", pressure_min);
            if (pressure_min < pressure_check_limit) {
                log::printf("\n WARNING: min pressure lower than min pressure limit: %g Pa\n",
                            pressure_min);
                if (exit_on_low_pressure_warning)
                    break;
            }

            // flush log output to file for security
            log::flush();

            logwriter.output_diagnostics(nstep,
                                         simulation_time,
                                         total_bytes,
                                         free_bytes,
                                         elapsed_time,
                                         time_left,
                                         mean_delta_per_step,
                                         end_time);
            log::printf("--\n");
        }
        else {
            printf("\n Time step number = %d/%d || Time = %f days. \n\t Elapsed %s || Left: %s || "
                   "Completion: %s. %s",
                   nstep,
                   nsmax,
                   simulation_time / 86400.,
                   duration_to_str(elapsed_time).c_str(),
                   duration_to_str(time_left).c_str(),
                   end_time_str.str().c_str(),
                   file_output ? "[saved output]" : "");
#ifdef STEP_TIMING_INFO
            printf("\n\t Step duration: %f s || Mean step duration %f s",
                   this_step_delta,
                   mean_delta_per_step);
#endif // STEP_TIMING_INFO

            // Run pressure check
            double pressure_min = gpu_min_on_device<1024>(X.pressure_d, X.point_num * X.nv);
            printf("\n min pressure : %g Pa\n", pressure_min);
            if (pressure_min < pressure_check_limit) {
                printf("\n WARNING: min pressure lower than min pressure limit: %g Pa\n",
                       pressure_min);
                if (exit_on_low_pressure_warning)
                    break;
            }
            printf("--\n");
        }


        if (caught_signal != ESIG_NOSIG) {
            //exit loop and application after save on SIGTERM or SIGINT
            break;
        }
        n_since_output++;
        // Marker to see to what step this output belongs
    }
    //
    //  Prints the duration of the integration.
    long finishTime = clock();
    log::printf("\n\n Integration time = %f seconds\n\n",
                double((finishTime - startTime)) / double(CLOCKS_PER_SEC));

    //
    //  Checks for errors in the device.
    cudaDeviceSynchronize();
    error = cudaGetLastError();
    log::printf("CudaMalloc error = %d = %s\n\n", error, cudaGetErrorString(error));
    //
    //  END OF THE ESP
    log::printf("End of the ESP!\n\n");

    cuda_device_memory_manager::get_instance().deallocate();

    Grid.free_memory();

    return 0;
}
