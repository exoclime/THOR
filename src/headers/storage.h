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
// Description: Store binary arrays to file
//
//
//
// Method: Write to HDF5 files
//
// Known limitations: None.
//
//
// Known issues: None.
//
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

#include <iostream>
#include <memory>
#include <string>

#include "H5Cpp.h"

#include "log_writer.h"

using std::string;
using namespace H5;

#if H5_VERSION_GE(1,10,2)
#define PRINT_EXCEPTION(exception) exception.printErrorStack()
#else
#define PRINT_EXCEPTION(exception) exception.printError()
#endif


class storage
{
public:
    storage(const string& filename, const bool& read = false);


    bool has_table(string name) {
        if (file != nullptr) {
            try { // to determine if the dataset exists in the group
                DataSet dataset = file->openDataSet(name);
                return true;
            } catch (DataSetIException & error) {
                return false;
            } catch (FileIException & error) {
                return false;
            }
        }
        else
            return false;
    }


    // Store a table of type T - double or int, in output file
    template<typename T>
    void append_table(T* data, const int& size, string name, string unit, string description);

    // read table of type T from output file
    template<typename T> bool read_table(const string& name, std::unique_ptr<T[]>& data, int& size);

    template<typename T> bool read_table_to_ptr(const string& name, T* data, int size);


    // write a scalar value to output
    template<typename T> void append_value(T value, string name, string unit, string description) {
        append_table(&value, 1, name, unit, description);
    }

    // read a scalar value from input
    template<typename T> bool read_value(const string& name, T& data) {

        std::unique_ptr<T[]> buf  = nullptr;
        int                  size = 0;

        bool out = read_table(name, buf, size);

        if (out && size == 1)
            data = buf[0];
        if (size != 1)
            return false;

        return out;
    }


    // Template functions for data type detection in append_table function.
    template<typename T> DataType get_datatype(T& input) {
        log::printf("data type not supported for storage.\n");

        throw std::runtime_error("data type not supported for storage");

        //       return PredType::STD_REF_OBJ;
    }

    DataType get_datatype(double& input) {
        return PredType::IEEE_F64LE;
    }

    DataType get_datatype(int& input) {
        return PredType::STD_I32LE;
    }

private:
    std::unique_ptr<H5File> file;
};

template<typename T>
void storage::append_table(T* data, const int& size, string name, string unit, string description) {
    if (file != nullptr) {
        try {

            DataType dt = get_datatype(data[0]);


            //create dataspace
            hsize_t   dims[] = {(hsize_t)size}; // dimensions of dataset
            DataSpace dataspace(1, dims);

            // create dataset with default properties
            DataSet dataset(file->createDataSet("/" + name, dt, dataspace));

            dataset.write(data, dt);

            {
                StrType   strdatatype(PredType::C_S1, description.length());
                DataSpace attrspace = H5::DataSpace(H5S_SCALAR);
                Attribute attr      = dataset.createAttribute("Variable", strdatatype, attrspace);
                attr.write(strdatatype, description.c_str());
            }
            {
                StrType   strdatatype(PredType::C_S1, unit.length());
                DataSpace attrspace = H5::DataSpace(H5S_SCALAR);
                Attribute attr      = dataset.createAttribute("units", strdatatype, attrspace);
                attr.write(strdatatype, unit.c_str());
            }
        } // end of try block
        // catch failure caused by the H5File operations
        catch (FileIException & error) {
            log::printf("FileIException: %s\t%d\t%d\n", name.c_str(), data, size);

            PRINT_EXCEPTION(error);
        }
        // catch failure caused by the DataSet operations
        catch (DataSetIException & error) {
            log::printf("DataSetIException: %s\t%d\t%d\n", name.c_str(), data, size);

            PRINT_EXCEPTION(error);
        }
        // catch failure caused by the DataSpace operations
        catch (DataSpaceIException & error) {
            log::printf("DataSpaceIException: %s\t%d\t%d\n", name.c_str(), data, size);

            PRINT_EXCEPTION(error);
        }
        // catch failure caused by the DataSpace operations
        catch (DataTypeIException & error) {
            log::printf("DataTypeIException: %s\t%d\t%d\n", name.c_str(), data, size);

            PRINT_EXCEPTION(error);
        }
    }
}

template<typename T>
bool storage::read_table(const string& name, std::unique_ptr<T[]>& data, int& size) {
    size = 0;

    if (file != nullptr) {
        DataType dt = get_datatype(data[0]);

        try {

            DataSet dataset = file->openDataSet(name);

            // get type in dataset
            DataType type = dataset.getDataType();


            if (type == dt) {
                // log::printf( "Data set has correct type.\n");
            }
            else {
                data = nullptr;
                size = 0;

                return false;
            }

            // get dataspace
            DataSpace dataspace = dataset.getSpace();
            // get dimensions and rank
            int     rank        = dataspace.getSimpleExtentNdims();
            hsize_t dims_out[1] = {0};
            if (rank == 1) {


                dataspace.getSimpleExtentDims(dims_out, NULL);
            }
            else {
                data = nullptr;
                size = 0;

                return false;
            }
            size = int(dims_out[0]);

            // build output array
            data = std::unique_ptr<T[]>(new T[(int)dims_out[0]], std::default_delete<T[]>());


            dataset.read(data.get(), dt);

        }

        // end of try block
        // catch failure caused by the H5File operations
        catch (FileIException & error) {
            PRINT_EXCEPTION(error);
            data = nullptr;

            return false;
        }
        // catch failure caused by the DataSet operations
        catch (DataSetIException & error) {
            PRINT_EXCEPTION(error);
            data = nullptr;

            return false;
        }
        // catch failure caused by the DataSpace operations
        catch (DataSpaceIException & error) {
            PRINT_EXCEPTION(error);
            data = nullptr;

            return false;
        }
        // catch failure caused by the DataSpace operations
        catch (DataTypeIException & error) {
            PRINT_EXCEPTION(error);
            data = nullptr;
            size = 0;

            return false;
        }
    }
    else {
        data = nullptr;
        size = 0;

        return false;
    }

    return true;
}

template<typename T> bool storage::read_table_to_ptr(const string& name, T* data, int size) {

    if (file != nullptr) {
        DataType dt = get_datatype(data[0]);

        try {

            DataSet dataset = file->openDataSet(name);

            // get type in dataset
            DataType type = dataset.getDataType();


            if (type == dt) {
                // log::printf( "Data set has correct type.\n");
            }
            else {
                return false;
            }

            // get dataspace
            DataSpace dataspace = dataset.getSpace();
            // get dimensions and rank
            int     rank        = dataspace.getSimpleExtentNdims();
            hsize_t dims_out[1] = {0};
            if (rank == 1) {
                dataspace.getSimpleExtentDims(dims_out, NULL);
            }
            else {
                return false;
            }

            if (size == int(dims_out[0])) {
                // build output array


                dataset.read(data, dt);
            }
            else
                return false;

        }

        // end of try block
        // catch failure caused by the H5File operations
        catch (FileIException & error) {
            PRINT_EXCEPTION(error);

            return false;
        }
        // catch failure caused by the DataSet operations
        catch (DataSetIException & error) {
            PRINT_EXCEPTION(error);

            return false;
        }
        // catch failure caused by the DataSpace operations
        catch (DataSpaceIException & error) {
            PRINT_EXCEPTION(error);

            return false;
        }
        // catch failure caused by the DataSpace operations
        catch (DataTypeIException & error) {
            PRINT_EXCEPTION(error);

            return false;
        }
    }
    else {

        return false;
    }

    return true;
}
