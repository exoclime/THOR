

// datatype to store values to update on device in THOR loop
#pragma once
#include <vector>

struct device_RK_array {
    double* array_d;
    double* arrayk_d;
    double* arrayi_d;
    int     dimensions;
    device_RK_array(){};

    device_RK_array(double* array_d_,
                    double* arrayk_d_,
                    double* arrayi_d_,
                    int     dimensions_) :
        array_d(array_d_),
        arrayk_d(arrayk_d_),
        arrayi_d(arrayi_d_),
        dimensions(dimensions_){};
};

class device_RK_array_manager
{
public:
    device_RK_array_manager();

    bool register_array(double* array_d, double* arrayk_d, double* arrayi, int dimensions);

    void allocate_device_array();

private:
    std::vector<device_RK_array> data;
};
