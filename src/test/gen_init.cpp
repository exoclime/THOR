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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "grid.h"
#include "planet.h"
#include <string>
#include <iostream>

#include "hdf5.h"
#include "storage.h"
#include <cstring>
#include <cmath>

using namespace std;

struct vector3
{
    double x;
    double y;
    double z;

    friend vector3 operator/(const vector3 & v, const double & scalar);
    friend vector3 operator*(const double & scalar,  const vector3 & v);
    friend double operator*(const vector3 & A, const vector3 & B);
    friend vector3 operator+(const vector3 & A, const vector3 & B);
    friend vector3 operator-(const vector3 & A, const vector3 & B);

    double norm() 
    {
        return sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
    }    
};


vector3 operator/(const vector3 & v, const double & scalar)
{
    return vector3({v.x/scalar, v.y/scalar, v.z/scalar});
}

vector3 operator*(const double & scalar,  const vector3 & v)
{
    return vector3({v.x*scalar, v.y*scalar, v.z*scalar});
}
double operator*(const vector3 & A, const vector3 & B)
{
    return A.x*B.x + A.y*B.y + A.z*B.z;
}

    
vector3 operator+(const vector3 & A, const vector3 & B)
{
    return vector3({A.x+B.x,A.y+B.y,A.z+B.z});
}

            
vector3 operator-(const vector3 & A, const vector3 & B)
{
    return vector3({A.x-B.x,A.y-B.y,A.z-B.z});
}

struct basis_sph
{
    vector3 e_r;
    vector3 e_theta;
    vector3 e_phi;
};

basis_sph base_sph_vectors(double theta, double phi)
{
    
    basis_sph basis;
    basis.e_r.x = cos(theta)*cos(phi);
    basis.e_r.y = cos(theta)*sin(phi);
    basis.e_r.z = sin(theta);

    basis.e_theta.x = -sin(theta)*cos(phi);
    basis.e_theta.y = -sin(theta)*sin(phi);
    basis.e_theta.z =  cos(theta);
    
    basis.e_phi.x = -sin(phi);
    basis.e_phi.y = cos(phi);
    basis.e_phi.z = 0.0;
    
   
   
    
    return basis;
}

void Output(int    ntstep         , // Number of integration steps
            int    fidx           , // Index of output file
            double Cp             , // Specific heat capacities [J/(Kg K)]
            double Rd             , // Gas constant [J/(Kg K)]
            double Omega          , // Rotation rate [s-1]
            double Gravit         , // Gravitational acceleration [m/s2]
            double Mmol           , // Mean molecular mass of dry air [kg]
            double P_Ref          , // Reference surface pressure [Pa] 
            double Top_altitude   , // Top of the model's domain [m]
            double A              , // Planet radius [m]
            bool spring_dynamics  ,
            double spring_beta    ,
            string simulation_ID  , // Planet ID
            double simulation_time, // Option for deep atmosphere
            const std::string & output_dir,
            XPlanet & planet,
            Icogrid & Grid,
            double * Rho_h,
            double * pressure_h,
            double * Mh_h,
            double * Wh_h){

//
//  Description: Model output.
//
    hid_t       file_id, dataset_id, att, dataspace_id;
    hid_t       stringType, stringSpace;
    hsize_t     dims[1];

    char FILE_NAME1[160];
    
    stringType =  H5Tcopy(H5T_C_S1);
    stringSpace=  H5Screate(H5S_SCALAR);
          
//  GRID OUTPUT
    if(ntstep == 0){           
        sprintf(FILE_NAME1, "%s/esp_output_grid_%s.h5", output_dir.c_str(), simulation_ID.c_str());

//      Create a new file using default properties. 
        file_id = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
//      Create the data space for the dataset. 
        dims[0] = Grid.nv; 
        dataspace_id = H5Screate_simple(1, dims, NULL);
//      Create the dataset. 
        dataset_id = H5Dcreate2(file_id, "/Altitude", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      Write the dataset. 
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Grid.Altitude);
//      Write attributes
        H5Tset_size(stringType, strlen("Altitude"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Altitude");
        H5Aclose(att);
        H5Tset_size(stringType, strlen("m"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "m");
//      End access to the dataset and release resources used by it. 
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);

//      Altitudeh
        dims[0] = Grid.nv+1; 
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/Altitudeh", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Grid.Altitudeh);
        H5Tset_size(stringType, strlen("Altitude at the interfaces"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Altitude at the interfaces");
        H5Aclose(att);
        H5Tset_size(stringType, strlen("m"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "m");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);

//      AreasT  
        dims[0] = Grid.point_num; 
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/areasT", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Grid.areasT);    
        H5Tset_size(stringType, strlen("Main cells areas"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Main cells areas");
        H5Aclose(att);
        H5Tset_size(stringType, strlen("m^2"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "m^2");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);         

//      Lon-lat grid
        dims[0] = 2*Grid.point_num;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/lonlat", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Grid.lonlat);    
        H5Tset_size(stringType, strlen("Longitudes and latitudes"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Longitudes and latitudes");
        H5Aclose(att);
        H5Tset_size(stringType, strlen("-"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "-");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);

//      Number of horizontal points
        double point_num_a[] = {(double)Grid.point_num};
        dims[0] = 1;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/point_num", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, point_num_a);    
        H5Tset_size(stringType, strlen("Number of grid points in one level"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Number of grid points in one level");
        H5Aclose(att);        
        H5Tset_size(stringType, strlen("-"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "-");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);

//      Number of vertical layers
        double nv_a[] = {(double)Grid.nv};
        dims[0] = 1;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/nv", H5T_IEEE_F64LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nv_a);    
        H5Tset_size(stringType, strlen("Number of vertical layers"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Number of vertical layers");
        H5Aclose(att);        
        H5Tset_size(stringType, strlen("-"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "-");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);     

//      point neighbours
        dims[0] = 6*Grid.point_num;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/pntloc", H5T_STD_I32LE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Grid.point_local);    
        H5Tset_size(stringType, strlen("Neighbours indexes"));
        att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "Neighbours indexes");
        H5Aclose(att);
        H5Tset_size(stringType, strlen("-"));
        att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(att, stringType, "-");
        H5Dclose(dataset_id);
        H5Aclose(att);
        H5Sclose(dataspace_id);
        
        
//      Close the file.
        H5Fflush(file_id, H5F_SCOPE_LOCAL);
        H5Fclose(file_id);
    }
  
//  PLANET
    if(ntstep == 0){           
        sprintf(FILE_NAME1, "%s/esp_output_planet_%s.h5", output_dir.c_str(), simulation_ID.c_str());
        file_id = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // glevel
        write_double_value_to_h5file(file_id, "/glevel", Grid.grid_level, "Horizontal subdivision level", "-");
        // vlevel
        write_double_value_to_h5file(file_id, "/vlevel", Grid.nv, "Vertical subdivision level", "-");
        // spring_dynamics
        write_double_value_to_h5file(file_id, "/spring_dynamics", spring_dynamics?1.0:0.0, "Spring dynamics", "-");
        // spring beta
        write_double_value_to_h5file(file_id, "/spring_beta", spring_beta, "Spring Beta", "-");        
        //      A
        write_double_value_to_h5file(file_id, "/A", planet.A, "Planet radius", "m");        
        //      Rd
        write_double_value_to_h5file(file_id, "/Rd", planet.Rd, "Gas constant", "J/(Kg K)");
        //      Omega
        write_double_value_to_h5file(file_id, "/Omega", planet.Omega, "Rotation rate", "1/s");
        //      Gravit
        write_double_value_to_h5file(file_id, "/Gravit", planet.Gravit,"Surface gravity" , "m/s^2");  
        //      Mmol
        write_double_value_to_h5file(file_id, "/Mmol", planet.Mmol,"Mean molecular mass of dry air" , "kg");        
        //      P_Ref
        write_double_value_to_h5file(file_id, "/P_Ref", planet.P_Ref, "Reference pressure","Pa" );
        //      Top_altitude
        write_double_value_to_h5file(file_id, "/Top_altitude", planet.Top_altitude, "Top of the model's domain", "m");
        //      CP
        write_double_value_to_h5file(file_id, "/Cp", planet.Cp, "Specific heat capacity", "J/(Kg K)");  
            
        H5Fflush(file_id, H5F_SCOPE_LOCAL);
        H5Fclose(file_id);
    }       
       
//  ESP OUTPUT
    sprintf(FILE_NAME1, "%s/esp_output_%s_%d.h5", output_dir.c_str(), simulation_ID.c_str(), fidx);
    file_id = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // step index
    write_int_value_to_h5file(file_id, "/nstep", ntstep, "Step number", "-" );
    
//  Simulation time
    dims[0] = 1; 
    double simulation_time_a[]           = {simulation_time};
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/simulation_time", H5T_IEEE_F64LE, dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, simulation_time_a);    
    H5Tset_size(stringType, strlen("Simulation time"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Simulation time");
    H5Aclose(att);        
    H5Tset_size(stringType, strlen("s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);    
    
//  Rho
    dims[0] = Grid.nv*Grid.point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Rho", H5T_IEEE_F64LE, dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Rho_h);
    H5Tset_size(stringType, strlen("Density"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Density");
    H5Aclose(att);    
    H5Tset_size(stringType, strlen("kg/m^3"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg/m^3");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);           
       
//  Pressure
    dims[0] = Grid.nv*Grid.point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Pressure", H5T_IEEE_F64LE, dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure_h);
    H5Tset_size(stringType, strlen("Pressure"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Pressure");
    H5Aclose(att);    
    H5Tset_size(stringType, strlen("Pa"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Pa");    
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);                   
    
//  Mh
    dims[0] = Grid.nv*Grid.point_num*3;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Mh", H5T_IEEE_F64LE, dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Mh_h);    
    H5Tset_size(stringType, strlen("Horizontal Momentum"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Horizontal Momentum");
    H5Aclose(att);    
    H5Tset_size(stringType, strlen("kg m/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);        
          
//  Wh
    dims[0] = Grid.nvi*Grid.point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Wh", H5T_IEEE_F64LE, dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Wh_h);    
    H5Tset_size(stringType, strlen("Vertical Momentum"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Vertical Momentum");
    H5Aclose(att);    
    H5Tset_size(stringType, strlen("kg m/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id); 
    
//  Close the file.
    H5Fflush(file_id, H5F_SCOPE_LOCAL);
    H5Fclose(file_id);
}


    
int main ()
{
    {
        bool spring_dynamics = true;
        double spring_beta = 1.15;
        bool sponge = false;
        
        
        // level 3 grid
        Icogrid Grid(spring_dynamics , // Spring dynamics option
                     spring_beta     , // Parameter beta for spring dynamics 
                     4               , // Horizontal resolution level
                     32              , // Number of vertical layers
                     10              , 
                     6371000.0       , // Planet radius
                     36000.0         ,
                     sponge);// Top of the model's domain

        string output_dir = "./gen/results_res_4/";

        
         XPlanet Planet;
         
         
         double * Rho_h        = (double*)malloc(Grid.nv*Grid.point_num   * sizeof(double));
         double * pressure_h   = (double*)malloc(Grid.nv*Grid.point_num   * sizeof(double));
         double * temperature_h= (double*)malloc(Grid.nv*Grid.point_num   * sizeof(double));
         double * Mh_h         = (double*)malloc(Grid.nv*Grid.point_num*3 * sizeof(double));
         double * Wh_h         = (double*)malloc(Grid.nvi*Grid.point_num  * sizeof(double));

         // write some initial conditions
         for (int i = 0; i < Grid.point_num; i++ ){
             double lon = Grid.lonlat[2*i];
             double lat = Grid.lonlat[2*i+1];
             basis_sph local_basis = base_sph_vectors(lat, lon);
             double Ha = Planet.Rd * Planet.Tmean / Planet.Gravit;
             printf("e_r: %f %f %f\n", local_basis.e_r.x, local_basis.e_r.y, local_basis.e_r.z);
             printf("e_theta: %f %f %f\n", local_basis.e_theta.x, local_basis.e_theta.y, local_basis.e_theta.z);
             printf("e_phi: %f %f %f\n", local_basis.e_phi.x, local_basis.e_phi.y, local_basis.e_phi.z);
             for (int lev = 0; lev < Grid.nv; lev++ ){
                 double p = lat/(M_PI/2.0);
                 pressure_h[i*Grid.nv + lev] = p;
                 
                 
                 
                 // pressure_h[i*Grid.nv + lev] = Planet.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                 //temperature_h[i*Grid.nv + lev] = Planet.Tmean;
             }
             
             for (int lev = 0; lev < Grid.nv; lev++ ){            
//              Density [kg/m3]
                 Rho_h[i*Grid.nv + lev] = lon/(M_PI);
                 
//              Momentum [kg/m3 m/s]
                 double v_theta = 1.0;
                 double v_phi = 0.0;
                 double v_r = 0.0;
                 
                 vector3 v_h = v_theta*local_basis.e_theta+
                     v_phi*local_basis.e_phi;
                 printf("%f %f %f\n", v_h.x, v_h.y, v_h.z);
                 
                 
                 //vector3 v_v = v_r*local_basis.e_r;
                 // Horizontal momentum in cartesian coordinate
                 Mh_h[i*3*Grid.nv + 3*lev + 0] = v_h.x;
                 Mh_h[i*3*Grid.nv + 3*lev + 1] = v_h.y;
                 Mh_h[i*3*Grid.nv + 3*lev + 2] = v_h.z;
                 
//              Vertical momentum [kg/m3 m/s]         
                 //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                 Wh_h[i*(Grid.nv+1) + lev] = v_r;// Layers interface.
             }
             Wh_h[i*(Grid.nv + 1) + Grid.nv] = 0.0;
         }
         //
//  Writes initial conditions
    double simulation_time = 0.0;
    Output(0                   ,
           0   ,
           Planet.Cp           , // Specific heat capacity [J/(Kg K)]
           Planet.Rd           , // Gas constant [J/(Kg K)]
           Planet.Omega        , // Rotation rate [s-1]
           Planet.Gravit       , // Gravitational acceleration [m/s2]
           Planet.Mmol         , // Mean molecular mass of dry air [kg]
           Planet.P_Ref        , // Reference surface pressure [Pa]
           Planet.Top_altitude , // Top of the model's domain [m]
           Planet.A            , // Planet Radius [m]
           spring_dynamics     ,
           spring_beta         ,
           Planet.simulation_ID, // Simulation ID (e.g., "Earth")
           simulation_time,
           output_dir,
           Planet,
           Grid,
           Rho_h,
           pressure_h,
           Mh_h,
           Wh_h);// Time of the simulation [s]


            // write some initial conditions
         for (int i = 0; i < Grid.point_num; i++ ){
             double lon = Grid.lonlat[2*i];
             double lat = Grid.lonlat[2*i+1];
             basis_sph local_basis = base_sph_vectors(lat, lon);
             double Ha = Planet.Rd * Planet.Tmean / Planet.Gravit;
             for (int lev = 0; lev < Grid.nv; lev++ ){
                 double p = lat/(M_PI/2.0);
                 pressure_h[i*Grid.nv + lev] = p;
                 // pressure_h[i*Grid.nv + lev] = Planet.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                 // temperature_h[i*Grid.nv + lev] = Planet.Tmean;
             }
             
             for (int lev = 0; lev < Grid.nv; lev++ ){            
//              Density [kg/m3]
                 // Rho_h[i*Grid.nv + lev] = pressure_h[i*Grid.nv + lev] / (temperature_h[i*Grid.nv + lev] * Planet.Rd);
                 Rho_h[i*Grid.nv + lev] = lon/(M_PI);
                 
//              Momentum [kg/m3 m/s]
                 double v_theta = cos(lev/double(Grid.nv)*2.0*M_PI);
                 double v_phi = sin(lev/double(Grid.nv)*2.0*M_PI);
                 double v_r = 0.0;

                 
                 
                 vector3 v_h = v_theta*local_basis.e_theta+
                     v_phi*local_basis.e_phi ;
                 
                 //vector3 v_v = v_r*local_basis.e_r;
                 
                 Mh_h[i*3*Grid.nv + 3*lev + 0] = v_h.x;
                 Mh_h[i*3*Grid.nv + 3*lev + 1] = v_h.y;
                 Mh_h[i*3*Grid.nv + 3*lev + 2] = v_h.z;
                 
//              Vertical momentum [kg/m3 m/s]         
                 //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                 Wh_h[i*(Grid.nv+1) + lev] = v_r;// Layers interface.
             }
             Wh_h[i*(Grid.nv + 1) + Grid.nv] = 0.0;
         }

         //  Writes initial conditions
     simulation_time = 1.0;
    Output(1                   ,
           1   ,
           Planet.Cp           , // Specific heat capacity [J/(Kg K)]
           Planet.Rd           , // Gas constant [J/(Kg K)]
           Planet.Omega        , // Rotation rate [s-1]
           Planet.Gravit       , // Gravitational acceleration [m/s2]
           Planet.Mmol         , // Mean molecular mass of dry air [kg]
           Planet.P_Ref        , // Reference surface pressure [Pa]
           Planet.Top_altitude , // Top of the model's domain [m]
           Planet.A            , // Planet Radius [m]
           spring_dynamics     ,
           spring_beta         ,
           Planet.simulation_ID, // Simulation ID (e.g., "Earth")
           simulation_time,
           output_dir,
           Planet,
           Grid,
           Rho_h,
           pressure_h,
           Mh_h,
           Wh_h);// Time of the simulation [s]
            // write some initial conditions
         for (int i = 0; i < Grid.point_num; i++ ){
             double lon = Grid.lonlat[2*i];
             double lat = Grid.lonlat[2*i+1];
             basis_sph local_basis = base_sph_vectors(lat, lon);
             double Ha = Planet.Rd * Planet.Tmean / Planet.Gravit;
             for (int lev = 0; lev < Grid.nv; lev++ ){
                 double p = lat/(M_PI/2.0);
                 pressure_h[i*Grid.nv + lev] = p;
                 // pressure_h[i*Grid.nv + lev] = Planet.P_Ref*exp(-Grid.Altitude[lev] / Ha);
                 // temperature_h[i*Grid.nv + lev] = Planet.Tmean;
             }
             
             for (int lev = 0; lev < Grid.nv; lev++ ){            
//              Density [kg/m3]
                 Rho_h[i*Grid.nv + lev] = lon/(M_PI);
                 
                 // Rho_h[i*Grid.nv + lev] = pressure_h[i*Grid.nv + lev] / (temperature_h[i*Grid.nv + lev] * Planet.Rd);
                 
//              Momentum [kg/m3 m/s]
                 double phase = lev/double(Grid.nv)*2.0*M_PI;
                 
                 double v_theta = cos(lat + phase);
                 double v_phi = sin(lat + phase);
                 double v_r = 0.5;
                 

                 
                 
                 vector3 v_h = v_theta*local_basis.e_theta+
                     v_phi*local_basis.e_phi;
                 
                 //vector3 v_v = v_r*local_basis.e_r;
                 
                 Mh_h[i*3*Grid.nv + 3*lev + 0] = v_h.x;
                 Mh_h[i*3*Grid.nv + 3*lev + 1] = v_h.y;
                 Mh_h[i*3*Grid.nv + 3*lev + 2] = v_h.z;
                 
//              Vertical momentum [kg/m3 m/s]         
                 //W_h[i*Grid.nv + lev] = 0.0;     // Center of the layer.
                 Wh_h[i*(Grid.nv+1) + lev] = v_r;// Layers interface.
             }
             Wh_h[i*(Grid.nv + 1) + Grid.nv] = 0.0;
         }

         //  Writes initial conditions
    simulation_time = 2.0;
    Output(2                   ,
           2   ,
           Planet.Cp           , // Specific heat capacity [J/(Kg K)]
           Planet.Rd           , // Gas constant [J/(Kg K)]
           Planet.Omega        , // Rotation rate [s-1]
           Planet.Gravit       , // Gravitational acceleration [m/s2]
           Planet.Mmol         , // Mean molecular mass of dry air [kg]
           Planet.P_Ref        , // Reference surface pressure [Pa]
           Planet.Top_altitude , // Top of the model's domain [m]
           Planet.A            , // Planet Radius [m]
           spring_dynamics     ,
           spring_beta         ,
           Planet.simulation_ID, // Simulation ID (e.g., "Earth")
           simulation_time,
           output_dir,
           Planet,
           Grid,
           Rho_h,
           pressure_h,
           Mh_h,
           Wh_h);// Time of the simulation [s]

    
    }
    
    
    exit(0);
}
