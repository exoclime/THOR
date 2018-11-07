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
// Timer helper to measure tome for iterations, elapsed time and time remaining
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
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

#include "iteration_timer.h"


iteration_timer::iteration_timer(int initial_num_steps_, int max_steps_):
    max_steps(max_steps_),
    initial_num_steps(initial_num_steps_)

{
    start_sim = std::chrono::system_clock::now();
}

void iteration_timer::iteration(int          nstep,
                                double&      mean_delta_per_step,
                                double&      elapsed_time,
                                double&      time_left_,
                                std::time_t& end_time) {
    // Get current time
    std::chrono::system_clock::time_point end_step = std::chrono::system_clock::now();

    // time since simulation start
    std::chrono::duration<double, std::ratio<1L, 1L>> sim_delta = end_step - start_sim;

    // number of steps since simulation start and to end of simulation
    long num_steps_elapsed = nstep - initial_num_steps + 1;
    long num_steps_left    = max_steps - num_steps_elapsed;

    // mean length of step
    mean_delta_per_step = sim_delta.count() / double(num_steps_elapsed);

    // time left to end of simulation
    std::chrono::duration<double, std::ratio<1L, 1L>> time_left(double(num_steps_left) * mean_delta_per_step);

    // estimated time of simulation
    std::chrono::system_clock::time_point sim_end = end_step + std::chrono::duration_cast<std::chrono::microseconds>(time_left);
    // format output
    end_time = std::chrono::system_clock::to_time_t(sim_end);

    elapsed_time = sim_delta.count();
    time_left_   = time_left.count();
}
