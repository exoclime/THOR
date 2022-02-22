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
// Timer helper to measure time for iterations, elapsed time and time remaining
//

//
// Known limitations:
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

#pragma once

#include <chrono>

class iteration_timer
{
public:
    iteration_timer(int initial_num_steps, int max_steps);

    void iteration(int          nstep,
                   double&      mean_delta_per_step,
                   double&      step_delta,
                   double&      elapsed_time,
                   double&      time_left,
                   std::time_t& end_time);

private:
    int                                   max_steps         = 1;
    int                                   initial_num_steps = 0;
    std::chrono::system_clock::time_point start_sim;
    std::chrono::system_clock::time_point previous_sim_step;
};
