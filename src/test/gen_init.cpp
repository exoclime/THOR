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
// Description: generate grid and initial conditions
//
//
//
// Method: [1] - Dumps output to binary file on a flag
//         [2] - Reads data from binary files on a flag and compare to
//               dumped output
//
//
//
// Known limitations: None.
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

#include "grid.h"
#include "simulation_setup.h"
#include <iostream>
#include <string>

#include "hdf5.h"
#include "storage.h"
#include "directories.h"
#include <cmath>
#include <cstring>

using namespace std;

struct vector3 {
    double x;
    double y;
    double z;

    friend vector3 operator/(const vector3& v, const double& scalar);
    friend vector3 operator*(const double& scalar, const vector3& v);
    friend double  operator*(const vector3& A, const vector3& B);
    friend vector3 operator+(const vector3& A, const vector3& B);
    friend vector3 operator-(const vector3& A, const vector3& B);

    double norm() {
        return sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
    }
};


vector3 operator/(const vector3& v, const double& scalar) {
    return vector3({v.x / scalar, v.y / scalar, v.z / scalar});
}

vector3 operator*(const double& scalar, const vector3& v) {
    return vector3({v.x * scalar, v.y * scalar, v.z * scalar});
}
double operator*(const vector3& A, const vector3& B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
}


vector3 operator+(const vector3& A, const vector3& B) {
    return vector3({A.x + B.x, A.y + B.y, A.z + B.z});
}


vector3 operator-(const vector3& A, const vector3& B) {
    return vector3({A.x - B.x, A.y - B.y, A.z - B.z});
}

struct basis_sph {
    vector3 e_r;
    vector3 e_theta;
    vector3 e_phi;
};

basis_sph base_sph_vectors(double theta, double phi) {

    basis_sph basis;
    basis.e_r.x = cos(theta) * cos(phi);
    basis.e_r.y = cos(theta) * sin(phi);
    basis.e_r.z = sin(theta);

    basis.e_theta.x = -sin(theta) * cos(phi);
    basis.e_theta.y = -sin(theta) * sin(phi);
    basis.e_theta.z = cos(theta);

    basis.e_phi.x = -sin(phi);
    basis.e_phi.y = cos(phi);
    basis.e_phi.z = 0.0;


    return basis;
}


void Output(int                current_step, // Number of integration steps
            int                fidx,         // Index of output file
            bool               spring_dynamics,
            double             spring_beta,
            string             simulation_ID,   // Planet ID
            double             simulation_time, // Option for deep atmosphere
            const std::string& output_dir,
            SimulationSetup&   sim,
            Icogrid&           Grid,
            double*            Rho_h,
            double*            pressure_h,
            double*            Mh_h,
            double*            Wh_h) {

    //
    //  Description: Model output.
    //
    char FILE_NAME1[160];

    //  GRID OUTPUT
    if (current_step == 0) {
        sprintf(FILE_NAME1, "%s/esp_output_grid_%s.h5", output_dir.c_str(), simulation_ID.c_str());

        storage s(FILE_NAME1);

        s.append_table(Grid.Altitude,
                       Grid.nv,
                       "/Altitude",
                       "m",
                       "Altitude");

        //      Altitudeh
        s.append_table(Grid.Altitudeh,
                       Grid.nv + 1,
                       "/Altitudeh",
                       "m",
                       "Altitude at the interfaces");

        //      AreasT
        s.append_table(Grid.areasT,
                       Grid.point_num,
                       "/areasT",
                       "m^2",
                       "Main cells areas");

        //      Lon-lat grid
        s.append_table(Grid.lonlat,
                       2 * Grid.point_num,
                       "/lonlat",
                       "-",
                       "Longitudes and latitudes");

        //      Number of horizontal points
        s.append_value((double)Grid.point_num,
                       "/point_num",
                       "-",
                       "Number of grid points in one level");

        //      Number of vertical layers
        s.append_value((double)Grid.nv,
                       "/nv",
                       "-",
                       "Number of vertical layers");
        //      point neighbours
        s.append_table(Grid.point_local,
                       6 * Grid.point_num,
                       "/pntloc",
                       "-",
                       "Neighbours indexes");
    }

   
    
    //  PLANET
    if (current_step == 0) {
        sprintf(FILE_NAME1, "%s/esp_output_planet_%s.h5", output_dir.c_str(), simulation_ID.c_str());
        storage s(FILE_NAME1);

        // glevel
        s.append_value(Grid.grid_level, "/glevel", "-", "Horizontal subdivision level");
        // vlevel
        s.append_value(Grid.nv, "/vlevel", "-", "Vertical subdivision level");
        // spring_dynamics
        s.append_value(spring_dynamics ? 1.0 : 0.0, "/spring_dynamics", "-", "Spring dynamics");
        // spring beta
        s.append_value(spring_beta, "/spring_beta", "-", "Spring Beta");
        //      A
        s.append_value(sim.A, "/A", "m", "Planet radius");
        //      Rd
        s.append_value(sim.Rd, "/Rd", "J/(Kg K)", "Gas constant");
        //      Omega
        s.append_value(sim.Omega, "/Omega", "1/s", "Rotation rate");
        //      Gravit
        s.append_value(sim.Gravit, "/Gravit", "m/s^2", "Surface gravity");
        //      P_Ref
        s.append_value(sim.P_Ref, "/P_Ref", "Pa", "Reference pressure");
        //      Top_altitude
        s.append_value(sim.Top_altitude, "/Top_altitude", "m", "Top of the model's domain");
        //      CP
        s.append_value(sim.Cp, "/Cp", "J/(Kg K)", "Specific heat capacity");
        //      SpongeLayer option
        s.append_value(sim.SpongeLayer ? 1.0 : 0.0, "/SpongeLayer", "-", "Using SpongeLayer?");
        //      DeepModel option
        s.append_value(sim.DeepModel ? 1.0 : 0.0, "/DeepModel", "-", "Using Deep Model");

        //      NonHydro option
        s.append_value(sim.NonHydro ? 1.0 : 0.0, "/NonHydro", "-", "Using Non Hydrostatic parameter");

        //      DivDampP option
        s.append_value(sim.DivDampP ? 1.0 : 0.0, "/DivDampP", "-", "Using Divergence-damping");

        //      HyDiff option
        s.append_value(sim.HyDiff ? 1.0 : 0.0, "/HyDiff", "-", "Using Hyper diffusion");

        //      Hyperdiffusion strength
        s.append_value(sim.Diffc, "/Diffc", "-", "Hyperdiffusion strength");

        //      Tmean
        s.append_value(sim.Tmean, "/Tmean", "K", "Mean atmospheric temperature");

        //      conv adj option
        s.append_value(sim.conv_adj, "/conv_adj", "-", "Using convection adjustment");

        //      GCM on/off option
        s.append_value(sim.gcm_off ? 1.0 : 0.0, "/gcm_off", "-", "GCM off");

        //      rest option
        s.append_value(sim.rest ? 1.0 : 0.0, "/rest", "-", "Starting from rest");

        //      TPprof option
        s.append_value(sim.TPprof, "/TPprof", "-", "Initial TP profile option");
    }

    //  ESP OUTPUT

    sprintf(FILE_NAME1, "%s/esp_output_%s_%d.h5", output_dir.c_str(), simulation_ID.c_str(), fidx);

    storage s(FILE_NAME1);
    // step index
    s.append_value(current_step,
                   "/nstep",
                   "-",
                   "Step number");

    //  Simulation time
    s.append_value(simulation_time,
                   "/simulation_time",
                   "s",
                   "Simulation time");

    //  Rho
    s.append_table(Rho_h,
                   Grid.nv * Grid.point_num,
                   "/Rho",
                   "kg/m^3",
                   "Density");

    //  Pressure
    s.append_table(pressure_h,
                   Grid.nv * Grid.point_num,
                   "/Pressure",
                   "Pa",
                   "Pressure");

    //  Mh
    s.append_table(Mh_h,
                   Grid.nv * Grid.point_num * 3,
                   "/Mh",
                   "kg m/s",
                   "Horizontal Momentum");

    //  Wh
    s.append_table(Wh_h,
                   Grid.nvi * Grid.point_num,
                   "/Wh",
                   "kg m/s",
                   "Vertical Momentum");
}


int main() {
    {
        bool   spring_dynamics = true;
        double spring_beta     = 1.15;
        bool   sponge          = false;


        int max_count;
        
        // level 3 grid
        Icogrid Grid(spring_dynamics, // Spring dynamics option
                     spring_beta,     // Parameter beta for spring dynamics
                     6,               // Horizontal resolution level
                     32,              // Number of vertical layers
                     10,
                     6371000.0, // Planet radius
                     36000.0,
                     sponge, // Top of the model's domain
                     &max_count);
        
        string output_dir = "./gen/results_res_6/";
        create_output_dir(output_dir);
        

        SimulationSetup sim;


        double* Rho_h         = (double*)malloc(Grid.nv * Grid.point_num * sizeof(double));
        double* pressure_h    = (double*)malloc(Grid.nv * Grid.point_num * sizeof(double));
        double* temperature_h = (double*)malloc(Grid.nv * Grid.point_num * sizeof(double));
        double* Mh_h          = (double*)malloc(Grid.nv * Grid.point_num * 3 * sizeof(double));
        double* Wh_h          = (double*)malloc(Grid.nvi * Grid.point_num * sizeof(double));

        // write some initial conditions
        for (int i = 0; i < Grid.point_num; i++) {
            double    lon         = Grid.lonlat[2 * i];
            double    lat         = Grid.lonlat[2 * i + 1];
            basis_sph local_basis = base_sph_vectors(lat, lon);
            double    Ha          = sim.Rd * sim.Tmean / sim.Gravit;
            printf("e_r: %f %f %f\n", local_basis.e_r.x, local_basis.e_r.y, local_basis.e_r.z);
            printf("e_theta: %f %f %f\n", local_basis.e_theta.x, local_basis.e_theta.y, local_basis.e_theta.z);
            printf("e_phi: %f %f %f\n", local_basis.e_phi.x, local_basis.e_phi.y, local_basis.e_phi.z);
            for (int lev = 0; lev < Grid.nv; lev++) {
                double p                      = lat / (M_PI / 2.0);
                pressure_h[i * Grid.nv + lev] = p;


                // pressure_h[i*Grid.nv + lev] = sim.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                //temperature_h[i*Grid.nv + lev] = sim.Tmean;
            }

            for (int lev = 0; lev < Grid.nv; lev++) {
                //              Density [kg/m3]
                Rho_h[i * Grid.nv + lev] = lon / (M_PI);

                //              Momentum [kg/m3 m/s]
                double v_theta = 1.0;
                double v_phi   = 0.0;
                double v_r     = 0.0;

                vector3 v_h = v_theta * local_basis.e_theta + v_phi * local_basis.e_phi;
                printf("%f %f %f\n", v_h.x, v_h.y, v_h.z);


                //vector3 v_v = v_r*local_basis.e_r;
                // Horizontal momentum in cartesian coordinate
                Mh_h[i * 3 * Grid.nv + 3 * lev + 0] = v_h.x;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 1] = v_h.y;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 2] = v_h.z;

                //              Vertical momentum [kg/m3 m/s]
                //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                Wh_h[i * (Grid.nv + 1) + lev] = v_r; // Layers interface.
            }
            Wh_h[i * (Grid.nv + 1) + Grid.nv] = 0.0;
        }
        //
        //  Writes initial conditions
        double simulation_time = 0.0;
        
        Output(0,
               0,
               spring_dynamics,
               spring_beta,
               sim.simulation_ID, // Simulation ID (e.g., "Earth")
               simulation_time,
               output_dir,
               sim,
               Grid,
               Rho_h,
               pressure_h,
               Mh_h,
               Wh_h); // Time of the simulation [s]


        // write some initial conditions
        for (int i = 0; i < Grid.point_num; i++) {
            double    lon         = Grid.lonlat[2 * i];
            double    lat         = Grid.lonlat[2 * i + 1];
            basis_sph local_basis = base_sph_vectors(lat, lon);
            double    Ha          = sim.Rd * sim.Tmean / sim.Gravit;
            for (int lev = 0; lev < Grid.nv; lev++) {
                double p                      = lat / (M_PI / 2.0);
                pressure_h[i * Grid.nv + lev] = p;
                // pressure_h[i*Grid.nv + lev] = sim.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                // temperature_h[i*Grid.nv + lev] = sim.Tmean;
            }

            for (int lev = 0; lev < Grid.nv; lev++) {
                //              Density [kg/m3]
                // Rho_h[i*Grid.nv + lev] = pressure_h[i*Grid.nv + lev] / (temperature_h[i*Grid.nv + lev] * sim.Rd);
                Rho_h[i * Grid.nv + lev] = lon / (M_PI);

                //              Momentum [kg/m3 m/s]
                double v_theta = cos(lev / double(Grid.nv) * 2.0 * M_PI);
                double v_phi   = sin(lev / double(Grid.nv) * 2.0 * M_PI);
                double v_r     = 0.0;


                vector3 v_h = v_theta * local_basis.e_theta + v_phi * local_basis.e_phi;

                //vector3 v_v = v_r*local_basis.e_r;

                Mh_h[i * 3 * Grid.nv + 3 * lev + 0] = v_h.x;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 1] = v_h.y;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 2] = v_h.z;

                //              Vertical momentum [kg/m3 m/s]
                //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                Wh_h[i * (Grid.nv + 1) + lev] = v_r; // Layers interface.
            }
            Wh_h[i * (Grid.nv + 1) + Grid.nv] = 0.0;
        }

        //  Writes initial conditions
        simulation_time = 1.0;
        Output(1,
               1,
               spring_dynamics,
               spring_beta,
               sim.simulation_ID, // Simulation ID (e.g., "Earth")
               simulation_time,
               output_dir,
               sim,
               Grid,
               Rho_h,
               pressure_h,
               Mh_h,
               Wh_h); // Time of the simulation [s]
                      // write some initial conditions
        for (int i = 0; i < Grid.point_num; i++) {
            double    lon         = Grid.lonlat[2 * i];
            double    lat         = Grid.lonlat[2 * i + 1];
            basis_sph local_basis = base_sph_vectors(lat, lon);
            double    Ha          = sim.Rd * sim.Tmean / sim.Gravit;
            for (int lev = 0; lev < Grid.nv; lev++) {
                double p                      = lat / (M_PI / 2.0);
                pressure_h[i * Grid.nv + lev] = p;
                // pressure_h[i*Grid.nv + lev] = sim.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                // temperature_h[i*Grid.nv + lev] = sim.Tmean;
            }

            for (int lev = 0; lev < Grid.nv; lev++) {
                //              Density [kg/m3]
                Rho_h[i * Grid.nv + lev] = lon / (M_PI);

                // Rho_h[i*Grid.nv + lev] = pressure_h[i*Grid.nv + lev] / (temperature_h[i*Grid.nv + lev] * sim.Rd);

                //              Momentum [kg/m3 m/s]
                double phase = lev / double(Grid.nv) * 2.0 * M_PI;

                double v_theta = cos(4.0 * lon) * cos(lat + phase);
                double v_phi   = cos(4.0 * lon) * sin(lat + phase);
                double v_r     = 0.5;


                vector3 v_h = v_theta * local_basis.e_theta + v_phi * local_basis.e_phi;

                //vector3 v_v = v_r*local_basis.e_r;

                Mh_h[i * 3 * Grid.nv + 3 * lev + 0] = v_h.x;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 1] = v_h.y;
                Mh_h[i * 3 * Grid.nv + 3 * lev + 2] = v_h.z;

                //              Vertical momentum [kg/m3 m/s]
                //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                Wh_h[i * (Grid.nv + 1) + lev] = v_r; // Layers interface.
            }
            Wh_h[i * (Grid.nv + 1) + Grid.nv] = 0.0;
        }

        //  Writes initial conditions
        simulation_time = 2.0;
        Output(2,
               2,
               spring_dynamics,
               spring_beta,
               sim.simulation_ID, // Simulation ID (e.g., "Earth")
               simulation_time,
               output_dir,
               sim,
               Grid,
               Rho_h,
               pressure_h,
               Mh_h,
               Wh_h); // Time of the simulation [s]
    }


    exit(0);
}
