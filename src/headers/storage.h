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

#pragma once

#include <string>
#include <memory>
#include <iostream>

#include "H5Cpp.h"

using std::cout;
using std::endl;

using std::string;
using namespace H5;


bool write_double_value_to_h5file(hid_t       file_id,
                                   const string & tablename,
                                  const double & out_value,
                                  const string & name,
                                  const string & unit );

    
bool write_double_table_to_h5file(hid_t       file_id,
                                  const string & tablename,
                                  double * double_table,
                                  int size,
                                  const string & name,
                                  const string & unit);

    
bool load_double_table_from_h5file(hid_t       file_id,
                                   const string & tablename,
                                   double * double_table,
                                   int expected_size );

bool load_double_value_from_h5file(hid_t       file_id,
                                   const string & tablename,
                                   double & out_value);

bool write_int_value_to_h5file(hid_t       file_id,
                               const string & tablename,
                               const int & out_value,
                               const string & name,
                               const string & unit );

    
bool write_int_table_to_h5file(hid_t       file_id,
                               const string & tablename,
                               int * double_table,
                               int size,
                               const string & name,
                               const string & unit);

    
bool load_int_table_from_h5file(hid_t       file_id,
                                const string & tablename,
                                int * double_table,
                                int expected_size );

bool load_int_value_from_h5file(hid_t       file_id,
                                const string & tablename,
                                int & out_value);
    
class storage
{
public:
    storage(const string & filename, const bool & read = false);
    

    // Store a table of type T - double or int, in output file
    template<typename T>
    void append_table(T * data,
                   const int & size,
                   string name,
                   string unit);

    // read table of type T from output file
    template<typename T>
    bool read_table(const string & name,
                    std::unique_ptr<T[]> & data,
                    int & size);

    // write a scalar value to output
    template<typename T>
    void append_value(T value,
                      string name,
                      string unit)
    {
        append_value(&value, 1, name, unit);
    }

    // read a scalar value from input
    template<typename T>
    void read_value(const string & name,
                    T & data)
    {
        std::unique_ptr<T[]> buf = &data;
        
        read_table(name, buf, 1);
    }
    
    
    // Template functions for data type detection in append_table function.
    template<typename T>
    DataType get_datatype(T & input)
    {
        cout << "data type not supported for storage" << endl;
        
        throw std::runtime_error("data type not supported for storage");
        
        //       return PredType::STD_REF_OBJ;
    }
    
    DataType get_datatype(double & input)
    {
        return PredType::IEEE_F64LE;
    }
    
    DataType get_datatype(int & input)
    {
        return PredType::STD_I32LE;
    }

private:
    std::unique_ptr<H5File> file;
    
                    
};

template<typename T>
void storage::append_table(T * data,
                   const int & size,
                   string name,
                   string unit)
{
    if (file != nullptr)
    {
        try
        {

            DataType dt = get_datatype(data[0]);
            
        
            //create dataspace
            hsize_t dims[] = {(hsize_t)size}; // dimensions of dataset
            DataSpace dataspace( 1, dims );
            
            // create dataset with default properties
            DataSet dataset(file->createDataSet("/" + name,
                                                dt,
                                                dataspace));

            dataset.write(data,  dt);

            {          
                StrType strdatatype(PredType::C_S1, name.length());
                hsize_t dims[] = {1}; // dimensions of attribute
                DataSpace attrspace = H5::DataSpace(1, dims);
                Attribute attr = dataset.createAttribute("Variable",
                                                         strdatatype,
                                                         attrspace);
                attr.write(strdatatype, name.c_str());
            }
            {          
                StrType strdatatype(PredType::C_S1, unit.length());
                hsize_t dims[] = {1}; // dimensions of attribute
                DataSpace attrspace = H5::DataSpace(1, dims);
                Attribute attr = dataset.createAttribute("units",
                                                         strdatatype,
                                                         attrspace);
                attr.write(strdatatype, unit.c_str());
            }            
        }  // end of try block
        // catch failure caused by the H5File operations
        catch( FileIException error )
        {
            cout << "FileIException: " << name
                 << "\t" << data
                 << "\t" << size << endl;
            
            error.printError();
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            cout << "DataSetIException: " << name
                 << "\t" << data
                 << "\t" << size << endl;

            error.printError();      
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            cout << "DataSpaceIException: " << name
                 << "\t" << data
                 << "\t" << size << endl;

            error.printError();
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            cout << "DataTypeIException: " << name
                 << "\t" << data
                 << "\t" << size << endl;

            error.printError();
        }
    }
}

template<typename T>
bool storage::read_table(const string & name,
                    std::unique_ptr<T[]> & data,
                    int & size)
{
    size = 0;
    
    if (file != nullptr)
    {
        DataType dt = get_datatype(data[0]);
        
        try
        {
        
            DataSet dataset = file->openDataSet( name );

            // get type in dataset
            DataType type = dataset.getDataType();

                  
            if( type == dt )
            {
                // cout << "Data set has correct type" << endl;
            }
            else
            {
                data = nullptr;
                size = 0;
                
                return false;
            }
            
            // get dataspace
            DataSpace dataspace = dataset.getSpace();
            // get dimensions and rank
            int rank = dataspace.getSimpleExtentNdims();
            hsize_t dims_out[1] = {0};
            if (rank == 1)
            {
                
                
                dataspace.getSimpleExtentDims( dims_out, NULL);
                // cout << "rank " << rank << ", dimensions " <<
                //    (unsigned long)(dims_out[0]) << " x " << endl;
            }
            else
            {
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
        catch( FileIException error )
        {
            error.printError();
            data = nullptr;
        
            return false;
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            error.printError();
            data = nullptr;
        
            return false;
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            error.printError();
            data = nullptr;
        
            return false;
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            error.printError();
            data = nullptr;
            size = 0;
            
            return false;
        }
    }
    else
    {
        data = nullptr;
        size = 0;
        
        return false;
    }

    return true;
    
    
}
