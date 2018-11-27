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


#include <iostream>

#include "storage.h"
#include <memory>

using namespace std;


const int num_d = 1024;
const int num_i = 100;


int main() {
    cout << "Storage test" << endl;

    {

        storage f("out.h5");

        uint32_t s = num_d;

        double d[s];

        for (int i = 0; i < num_d; i++)
            d[i] = double(i) / double(s);

        f.append_table(d, s, "Numbers", "m", "Number table");

        uint32_t s2 = num_i;

        int d2[s2];
        for (int i = 0; i < num_i; i++)
            d2[i] = i;

        f.append_table(d2, s2, "Indices", "m^2", "indices table");
    }

    {

        storage f("out.h5", true);

        {

            int size_out = 0;

            std::unique_ptr<double[]> data_ptr = nullptr;

            f.read_table("Numbers", data_ptr, size_out);
            if (size_out != num_d)
                cout << "error on size of table. Got: " << size_out
                     << "\texpected: " << num_d << endl;
            else {
                int cnt = 0;

                for (int i = 0; i < size_out; i++) {
                    double dat = double(i) / double(size_out);

                    if (data_ptr[i] != dat) {

                        cout << "wrong data(" << i << ")"
                             << "\tGot:\t"
                             << data_ptr[i] << "\texpected:\t" << dat << endl;
                        cnt++;
                    }
                }
                if (cnt > 0)
                    cout << "got " << cnt << " wrong data points" << endl;
            }
        }

        {
            int size_out = 0;

            std::unique_ptr<int[]> data_ptr = nullptr;

            f.read_table("Indices", data_ptr, size_out);
            if (size_out != num_i)
                cout << "error on size of table. Got: " << size_out
                     << "\texpected: " << num_i << endl;
            else {
                int cnt = 0;

                for (int i = 0; i < size_out; i++) {
                    if (data_ptr[i] != i)
                        cout << "wrong data(" << i << ")"
                             << "\tGot:\t"
                             << data_ptr[i] << "\texpected:\t" << i << endl;
                }
                if (cnt > 0)
                    cout << "got " << cnt << " wrong data points" << endl;
            }
        }
    }
}
