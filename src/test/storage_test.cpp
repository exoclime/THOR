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
// Description: test for storage class
//
//
// Method: writes an array to a file and reloads it
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



#include <iostream>

#include "storage.h"
#include <memory>

using namespace std;



int main()
{
    cout << "Storage test" << endl;

    {
        
        storage f("out.h5");

        uint32_t s = 1024;
    
        double d[s];
        
        for (int i = 0; i < s; i++)
            d[i] = double(i)/double(s);

        f.append_table(d, s, "Numbers", "m");
    }

    {
        
        storage f("out.h5" , true);

        int size_out = 0;

        std::unique_ptr<double[]> data_ptr = nullptr;
    
        f.read_table("Numbers", data_ptr, size_out);
        for (int i = 0; i < size_out; i++)
            cout << i << "\t"<< data_ptr[i]<<endl;
    }

}
