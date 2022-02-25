#pragma once

#include "cuda_device_memory.h"
#include "debug.h"
#include <string>

// forward declare ESP
class ESP;

using std::string;

struct diag_data {
    unsigned int flag;
    double4      data;
};


// Build a list of flags with a value for each bit, so that we can uses a binary mask
typedef enum {
    NAN_VALUE           = 1 << 0,
    NEGATIVE_VALUE      = 1 << 1,
    THOMAS_NOT_DD       = 1 << 2,
    THOMAS_BAD_SOLUTION = 1 << 3,
    BL_THOMAS_NOT_DD    = 1 << 4
} sim_error_flag;

const int DIAG_NUM_FLAGS = 3;

class kernel_diagnostics
{
public:
    kernel_diagnostics(ESP& esp);

    // The device memory for diagnostics
    // array of diag_data structure, contains diagnostics per element
    cuda_device_memory<diag_data> diagnostics;
    // globaal flag for the whole grid
    cuda_device_memory<unsigned int> diagnostics_global_flag;
    // set the global flag to 0, to use before calling a kernel to diagnose
    void reset_flag();
    // reset thhe global flag and the diagnostic memory for the full grid
    void reset();

    // check if a flag is enabled in global flag
    bool check_flag();
    // return global flag
    unsigned int get_flag();

    // dump data to file
    // checks the diagnostics array, element by element, to see if flags are set and print out flag and values.
    // Iterates on the full array.
    // specify the size of the grid, (as number of columns and number of levels)
    // that is used for diagnostics, depending if diagnosing a kernel computing one column value (point_num, 1),
    // or a volume element value (point_num, nv), or interface value (point_num, nvi)
    void dump_data(string name,              // base name of the checl, used to generate filename
                   int    iteration,         // iteration, for filename
                   int    subiteration,      // subiteration, for filename
                   int    subsubiteration,   // subsubiteration, for filename
                   ESP&   esp,               // esp reference
                   int    grid_num_points,   // size of grid
                   int    vertical_num_points); // number of levels
};
