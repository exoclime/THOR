#include "pylfrodull.h"
#include "alfrodullib.h"

#include <pybind11/pybind11.h>
#include <functional>

namespace py = pybind11;

py::function fn_clbck;
bool  fn_clbck_set = false;

void call_callback()
{
  printf("Called call_callback\n");
  fn_clbck();
  printf("Came back from callback\n");
  //  fn_clbck.release();
}

std::function<void()> calc_z_callback = call_callback;

void set_callback(py::object & fn)
{
  printf("Called set_callback\n");
  fn_clbck = fn;
  fn_clbck_set = true;
  
  set_z_calc_function(calc_z_callback);

}

auto cleanup_callback = []()
			{
			  printf("Cleanup pylfrodull module\n");
			  // perform cleanup here -- this function is called with the GIL held
			  if (fn_clbck_set)
			    {
			      fn_clbck.release();
			      fn_clbck_set = false;
			    }
};




PYBIND11_MODULE(pylfrodull, m) {


  
    m.doc() = "python wrapper for Alfrodull"; // optional module docstring

    m.def("pycompute_radiative_transfer", &wrap_compute_radiative_transfer, "compute radiative transfer");
    m.def("init_alfrodull", &init_alfrodull, "initialise Alfrodull Engine");
    m.def("deinit_alfrodull", &deinit_alfrodull, "deinitialise Alfrodull Engine");
    m.def("init_parameters", &init_parameters, "initialise global sim parameters");
    m.def("allocate", &allocate, "allocate internal memory");
    m.def("get_dev_pointers", &get_device_pointers_for_helios_write, "Get device pointers");
    m.def("prepare_planck_table", &prepare_planck_table, "Prepare planck table");
    m.def("correct_incident_energy", &correct_incident_energy, "Correct incident flux");

    m.def("get_opac_data_for_helios", &get_opac_data_for_helios, "get data for helios");
    m.def("set_clouds_data", &set_clouds_data, "set clouds data");
    // does not exist in this version anymore
    // m.def("set_surface_temperature", &set_surface_temperature, "set surface temperature");
   
    m.def("set_callback", &set_callback, "Set function callback");
    m.def("call_callback", &call_callback, "call function callback");
    m.add_object("_cleanup", py::capsule(cleanup_callback));
}
