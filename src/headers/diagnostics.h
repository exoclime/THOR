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
    THOMAS_BAD_SOLUTION = 1 << 3
} sim_error_flag;

const int DIAG_NUM_FLAGS = 3;

class kernel_diagnostics
{
public:
    kernel_diagnostics(ESP& esp);

    cuda_device_memory<diag_data>    diagnostics;
    cuda_device_memory<unsigned int> diagnostics_global_flag;
    void                             reset_flag();
    void                             reset();

    bool         check_flag();
    unsigned int get_flag();
    void         dump_data(string name,
                           int    iteration,
                           int    subiteration,
                           int    subsubiteration,
                           ESP&   esp,
                           int    grid_num_points,
                           int    vertical_num_points);
};
