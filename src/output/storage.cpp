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
// Description: Store binary arrays to a file, currently HDF5
//              abstracts storage type
//
//
//
// Method: Reads and write to files using names for output data
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


#include "storage.h"
#include "hdf5.h"
#include <iostream>

using namespace std;


storage::storage(const string& filename, const bool& read):
    file(nullptr) {
    if (read)
        file = std::unique_ptr<H5File>(new H5File(filename, H5F_ACC_RDONLY));
    else
        file = std::unique_ptr<H5File>(new H5File(filename, H5F_ACC_TRUNC));
}

// helpers of HDF5
// TODO: should end up being folded into storage class


// *******************************************************************
// * double table
// * write
bool write_double_table_to_h5file(hid_t         file_id,
                                  const string& tablename,
                                  double*       double_table,
                                  int           size,
                                  const string& name,
                                  const string& unit) {
    hid_t   dataset_id, att, dataspace_id;
    hid_t   stringType, stringSpace;
    hsize_t dims[1];

    stringType  = H5Tcopy(H5T_C_S1);
    stringSpace = H5Screate(H5S_SCALAR);

    dims[0]      = size;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id   = H5Dcreate2(file_id, tablename.c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, double_table);

    H5Tset_size(stringType, name.length());
    att = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, name.c_str());
    H5Aclose(att);
    H5Tset_size(stringType, unit.length());
    att = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, unit.c_str());
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

    return true;
}

bool write_double_value_to_h5file(hid_t         file_id,
                                  const string& tablename,
                                  const double& double_value,
                                  const string& name,
                                  const string& unit) {
    double val[] = {double_value};

    return write_double_table_to_h5file(file_id, tablename, val, 1, name, unit);
}


// * double table
// * read

bool load_double_table_from_h5file(hid_t         file_id,
                                   const string& tablename,
                                   double*       double_table,
                                   int           expected_size) {
    hid_t dataset_id;

    dataset_id = H5Dopen(file_id, tablename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {
        printf("Could not open dataset: %s.\n", tablename.c_str());
        return false;
    }

    hid_t   dspace     = H5Dget_space(dataset_id);
    hsize_t data_space = H5Sget_simple_extent_npoints(dspace);
    if (data_space != (unsigned int)(expected_size)) {
        printf("Initial condition load error for data: %s "
               "expected size: %u, got: %llu.\n",
               tablename.c_str(),
               (unsigned int)expected_size,
               data_space);

        return false;
    }
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, double_table);

    return true;
}

// * double table
// * write

bool load_double_value_from_h5file(hid_t         file_id,
                                   const string& tablename,
                                   double&       out_value) {
    return load_double_table_from_h5file(file_id, tablename, &out_value, 1);
}

// *******************************************************************
// * int table
// * write
bool write_int_table_to_h5file(hid_t         file_id,
                               const string& tablename,
                               int*          int_table,
                               int           size,
                               const string& name,
                               const string& unit) {
    hid_t   dataset_id, att, dataspace_id;
    hid_t   stringType, stringSpace;
    hsize_t dims[1];

    stringType  = H5Tcopy(H5T_C_S1);
    stringSpace = H5Screate(H5S_SCALAR);

    dims[0]      = size;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id   = H5Dcreate2(file_id, tablename.c_str(), H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, int_table);

    H5Tset_size(stringType, name.length());
    att = H5Acreate(dataset_id, "Variable", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, name.c_str());
    H5Aclose(att);
    H5Tset_size(stringType, unit.length());
    att = H5Acreate(dataset_id, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att, stringType, unit.c_str());
    H5Dclose(dataset_id);
    H5Aclose(att);
    H5Sclose(dataspace_id);

    return true;
}

bool write_int_value_to_h5file(hid_t         file_id,
                               const string& tablename,
                               const int&    int_value,
                               const string& name,
                               const string& unit) {
    int val[] = {int_value};

    return write_int_table_to_h5file(file_id, tablename, val, 1, name, unit);
}


// * int table
// * read

bool load_int_table_from_h5file(hid_t         file_id,
                                const string& tablename,
                                int*          int_table,
                                int           expected_size) {
    hid_t dataset_id;

    dataset_id = H5Dopen(file_id, tablename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {
        printf("Could not open dataset: %s.\n", tablename.c_str());
        return false;
    }

    hid_t   dspace     = H5Dget_space(dataset_id);
    hsize_t data_space = H5Sget_simple_extent_npoints(dspace);
    if (data_space != (unsigned int)(expected_size)) {
        printf("Initial condition load error for data: %s "
               "expected size: %u, got: %llu.\n",
               tablename.c_str(),
               (unsigned int)expected_size,
               data_space);

        return false;
    }
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, int_table);

    return true;
}

bool load_int_value_from_h5file(hid_t         file_id,
                                const string& tablename,
                                int&          out_value) {
    return load_int_table_from_h5file(file_id, tablename, &out_value, 1);
}
