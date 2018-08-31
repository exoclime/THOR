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
// Description: Writes the output.
//
//
// Method: Uses hdf5 files to write the output.
//
//
// Known limitations: None.
//
//
// Known issues: None.
//
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
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
#include <string>
#include "esp.h"
#include "hdf5.h"
#include "storage.h"




__host__ void ESP::CopyToHost(bool conservation){

//
//  Description: Transfer diagnostics from the device to the host.
//
        cudaMemcpy(Rho_h      , Rho_d      , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(Wh_h       , Wh_d       , point_num * nvi * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(pressure_h , pressure_d , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(Mh_h       , Mh_d       , 3 * point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
        if (conservation == true) {
          cudaMemcpy(Etotal_h   , Etotal_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(Mass_h     , Mass_d     , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(AngMomx_h  , AngMomx_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(AngMomy_h  , AngMomy_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(AngMomz_h  , AngMomz_d   , point_num * nv * sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(&GlobalE_h   , GlobalE_d    , sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(&GlobalMass_h, GlobalMass_d , sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(&GlobalAMx_h , GlobalAMx_d  , sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(&GlobalAMy_h , GlobalAMy_d  , sizeof(double), cudaMemcpyDeviceToHost);
          cudaMemcpy(&GlobalAMz_h , GlobalAMz_d  , sizeof(double), cudaMemcpyDeviceToHost);
        }
}

__host__ void ESP::Output(int    ntstep         , // Number of integration steps
                          int    fidx           , // Index of output file
                          double Cp             , // Specific heat capacities [J/(Kg K)]
                          double Rd             , // Gas constant [J/(Kg K)]
                          double Omega          , // Rotation rate [s-1]
                          double Gravit         , // Gravitational acceleration [m/s2]
                          double Mmol           , // Mean molecular mass of dry air [kg]
                          double P_Ref          , // Reference surface pressure [Pa]
                          double Top_altitude   , // Top of the model's domain [m]
                          double A              , // Planet radius [m]
                          char  *simulation_ID  , // Planet ID
                          double simulation_time, // Option for deep atmosphere
                          const std::string & output_dir,
                          bool conservation     ){

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
        sprintf(FILE_NAME1, "%s/esp_output_grid_%s.h5", output_dir.c_str(), simulation_ID);

//      Create a new file using default properties.
        file_id = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//      Create the data space for the dataset.
        dims[0] = nv;
        dataspace_id = H5Screate_simple(1, dims, NULL);
//      Create the dataset.
        dataset_id = H5Dcreate2(file_id, "/Altitude", H5T_IEEE_F64LE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      Write the dataset.
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Altitude_h);
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
        dims[0] = nv+1;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/Altitudeh", H5T_IEEE_F64LE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Altitudeh_h);
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
        dims[0] = point_num;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/areasT", H5T_IEEE_F64LE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, areasT_h);
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
        dims[0] = 2*point_num;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/lonlat", H5T_IEEE_F64LE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lonlat_h);
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
        double point_num_a[] = {(double)point_num};
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
        double nv_a[] = {(double)nv};
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
        dims[0] = 6*point_num;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate2(file_id, "/pntloc", H5T_STD_I32LE, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, point_local_h);
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
        sprintf(FILE_NAME1, "%s/esp_output_planet_%s.h5", output_dir.c_str(), simulation_ID);
        file_id = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // glevel
        write_double_value_to_h5file(file_id, "/glevel", glevel, "Horizontal subdivision level", "-");
        // vlevel
        write_double_value_to_h5file(file_id, "/vlevel", nv, "Vertical subdivision level", "-");
        // spring_dynamics
        write_double_value_to_h5file(file_id, "/spring_dynamics", spring_dynamics?1.0:0.0, "Spring dynamics", "-");
        // spring beta
        write_double_value_to_h5file(file_id, "/spring_beta", spring_beta, "Spring Beta", "-");
        //      A
        write_double_value_to_h5file(file_id, "/A", A, "Planet radius", "m");
        //      Rd
        write_double_value_to_h5file(file_id, "/Rd", Rd, "Gas constant", "J/(Kg K)");
        //      Omega
        write_double_value_to_h5file(file_id, "/Omega", Omega, "Rotation rate", "1/s");
        //      Gravit
        write_double_value_to_h5file(file_id, "/Gravit", Gravit,"Surface gravity" , "m/s^2");
        //      Mmol
        write_double_value_to_h5file(file_id, "/Mmol", Mmol,"Mean molecular mass of dry air" , "kg");
        //      P_Ref
        write_double_value_to_h5file(file_id, "/P_Ref", P_Ref, "Reference pressure","Pa" );
        //      Top_altitude
        write_double_value_to_h5file(file_id, "/Top_altitude", Top_altitude, "Top of the model's domain", "m");
        //      CP
        write_double_value_to_h5file(file_id, "/Cp", Cp, "Specific heat capacity", "J/(Kg K)");

        H5Fflush(file_id, H5F_SCOPE_LOCAL);
        H5Fclose(file_id);
    }

//  ESP OUTPUT
    sprintf(FILE_NAME1, "%s/esp_output_%s_%d.h5", output_dir.c_str(), simulation_ID, fidx);
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
    dims[0] = nv*point_num;
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
    dims[0] = nv*point_num;
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
    dims[0] = nv*point_num*3;
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
    dims[0] = nvi*point_num;
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

    if (conservation == true)  {
//  Etotal at each point
    dims[0] = nv*point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Etotal", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Etotal_h);
    H5Tset_size(stringType, strlen("Total Energy"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Total Energy");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s^2"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s^2");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  Mass at each point
    dims[0] = nv*point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Mass", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Mass_h);
    H5Tset_size(stringType, strlen("Mass"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Mass");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  AngMomx at each point
    dims[0] = nv*point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/AngMomx", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, AngMomx_h);
    H5Tset_size(stringType, strlen("AngMom in X"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "AngMom in X");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);
//
// //  AngMomy at each point
    dims[0] = nv*point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/AngMomy", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, AngMomy_h);
    H5Tset_size(stringType, strlen("AngMom in Y"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "AngMom in Y");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);
//
// //  AngMomz at each point
    dims[0] = nv*point_num;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/AngMomz", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, AngMomz_h);
    H5Tset_size(stringType, strlen("AngMom in Z"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "AngMom in Z");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  GlobalE (total energy over entire planet)
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/GlobalE", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GlobalE_h);
    H5Tset_size(stringType, strlen("Global Total Energy"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Global Total Energy");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s^2"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s^2");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  GlobalMass (total atmospheric mass over entire planet)
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/GlobalMass", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GlobalMass_h);
    H5Tset_size(stringType, strlen("Global Mass"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Global Mass");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  GlobalAMx (total angular momentum in x direction over entire planet)
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/GlobalAMx", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GlobalAMx_h);
    H5Tset_size(stringType, strlen("Global AngMomX"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Global AngMomX");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  GlobalAMy (total angular momentum in y direction over entire planet)
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/GlobalAMy", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GlobalAMy_h);
    H5Tset_size(stringType, strlen("Global AngMomY"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Global AngMomY");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

//  GlobalAMz (total angular momentum in y direction over entire planet)
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/GlobalAMz", H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GlobalAMz_h);
    H5Tset_size(stringType, strlen("Global AngMomZ"));
    att    = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "Global AngMomZ");
    H5Aclose(att);
    H5Tset_size(stringType, strlen("kg m^2/s"));
    att    = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, "kg m^2/s");
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);
    }

//  Close the file.
    H5Fflush(file_id, H5F_SCOPE_LOCAL);
    H5Fclose(file_id);
}
