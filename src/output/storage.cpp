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





storage::storage(const string & filename)
    : file(nullptr)
{
    file = std::make_unique<H5File>( filename, H5F_ACC_TRUNC );

    
}

storage::~storage()
{
    
}


void storage::append_table(double * data,
                   const int & size,
                   string name,
                   string unit)
{
    if (file != nullptr)
    {
        try
        {
        
            //create dataspace
            hsize_t dims[] = {(hsize_t)size}; // dimensions of dataset
            DataSpace dataspace( 1, dims );
            
            // create dataset with default properties
            DataSet dataset(file->createDataSet("/" + name,
                                                PredType::IEEE_F64LE,
                                                dataspace));

            dataset.write(data,  PredType::IEEE_F64LE);

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
            error.printError();
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            error.printError();      
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            error.printError();
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            error.printError();
        }
    }
}



void storage::write()
{

}


bool storage::read_table(string name,
                    float * data,
                    int & size)
{
    return false;
}
