#include "cuda_device_memory.h"


cuda_device_memory_manager& cuda_device_memory_manager::get_instance() {
  static cuda_device_memory_manager cdmm;
  
  return cdmm;
}

cuda_device_memory_manager::cuda_device_memory_manager() {
  
};

cuda_device_memory_manager::~cuda_device_memory_manager() {
  deallocate();
};

void cuda_device_memory_manager::register_mem(cuda_device_memory_interface * cdm)
{
  device_memory.push_back(cdm);
};

void cuda_device_memory_manager::unregister_mem(cuda_device_memory_interface * cdm)
{
  
  vector<cuda_device_memory_interface*>::iterator position = std::find(device_memory.begin(), device_memory.end(), cdm);
  if (position != device_memory.end()) // == myVector.end() means the element was not found
    device_memory.erase(position);
};

void cuda_device_memory_manager::deallocate()
{
  for (auto & cdm: device_memory)
    {
      cdm->deallocate();
    }
};

