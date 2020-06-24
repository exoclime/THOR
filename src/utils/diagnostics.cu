#include "diagnostics.h"
#include "directories.h"
#include "esp.h"
#include <memory>

kernel_diagnostics::kernel_diagnostics(ESP& esp) {
    // use number of points times number of interfaces, to fit all data sizes
    diagnostics.allocate(esp.point_num * esp.nvi);
    diagnostics_global_flag.allocate(1);
}

void kernel_diagnostics::reset_flag() {
    diagnostics_global_flag.zero();
}

void kernel_diagnostics::reset() {
    diagnostics_global_flag.zero();
    diagnostics.zero();
}

bool kernel_diagnostics::check_flag() {
    std::shared_ptr<unsigned int[]> global_flag_h = diagnostics_global_flag.get_host_data();
    return global_flag_h[0] != 0;
}

// Dump data for grid, following grid config passed as argument
void kernel_diagnostics::dump_data(string name,
                                   int    iteration,
                                   int    subiteration,
                                   int    subsubiteration,
                                   ESP&   esp,
                                   int    grid_num_points,
                                   int    vertical_num_points) {
    string output_dir = esp.get_output_dir();
    path   o(output_dir);
    o /= string("diagnostics");

    create_output_dir(o.to_string());
    o /= (string("diag_") + name + "_" + std::to_string(iteration) + "_"
          + std::to_string(subiteration) + "_" + std::to_string(subsubiteration) + ".txt");
    FILE* pFile = fopen(o.to_string().c_str(), "w");

    std::shared_ptr<diag_data[]> diag_data_h = diagnostics.get_host_data();

    for (int p = 0; p < grid_num_points; p++) {
        for (int lev = 0; lev < vertical_num_points; lev++) {
            int idx = p * vertical_num_points + lev;
            if (diag_data_h[idx].flag != 0) {
                if ((diag_data_h[idx].flag & NAN_VALUE) != 0)
                    fprintf(pFile, "NAN_VALUE - idx: %d - lev: %d\n", p, lev);
                if ((diag_data_h[idx].flag & NEGATIVE_VALUE) != 0)
                    fprintf(pFile,
                            "NEGATIVE_VALUE - idx: %d - lev: %d - value: %g\n",
                            p,
                            lev,
                            diag_data_h[idx].data.x);
                if ((diag_data_h[idx].flag & THOMAS_NOT_DD) != 0)
                    fprintf(pFile,
                            "THOMAS_NOT_DD - idx: %d - lev: %d - aa: %g - bb: %g - cc_s: %g - sum: "
                            "%g\n",
                            p,
                            lev,
                            diag_data_h[idx].data.x,
                            diag_data_h[idx].data.y,
                            diag_data_h[idx].data.z,
                            diag_data_h[idx].data.w);
                if ((diag_data_h[idx].flag & THOMAS_BAD_SOLUTION) != 0)
                    fprintf(pFile, "THOMAS_BAD_SOLUTION - idx: %d - lev: %d\n", p, lev);
            }
        }
    }
    fclose(pFile);
}
