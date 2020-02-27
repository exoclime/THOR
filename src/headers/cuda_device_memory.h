#pragma once

#include <cstddef>
#include <cstdio>
#include <memory>

#include <vector>

#include <algorithm>

using std::vector;

// class for memory manager storage, to be able to store multiple template
// instantiations
class cuda_device_memory_interface {
public:
  virtual void deallocate() = 0;
};


// singleton class for memory management
// necessary to free device memory before the whole cuda environment is tear down.
class cuda_device_memory_manager
{
public:
  
  // make a singleton, so that the object exists only once
  // use this to get a reference to the object
  static cuda_device_memory_manager& get_instance();

  ~cuda_device_memory_manager();
  
  // no copy constructor and assignement operator
  cuda_device_memory_manager(cuda_device_memory_manager const&) = delete;
  void operator=(cuda_device_memory_manager const&) = delete;

  void register_mem(cuda_device_memory_interface * cdm);

  void unregister_mem(cuda_device_memory_interface * cdm);
  

  void deallocate();

private:
  // make constructor private, can only be instantiated through get_instance
  cuda_device_memory_manager();
  
  vector<cuda_device_memory_interface*> device_memory;
};


// class to manage device memory, to take care of allocation and deallocation.
template<typename T> class cuda_device_memory : cuda_device_memory_interface
{
public:
    cuda_device_memory(){
      register_to_cdmm();
    };

    cuda_device_memory(size_t size) {
      register_to_cdmm();
      allocate(size);
    };

    ~cuda_device_memory() {
      deallocate();
      unregister_from_cdmm();
    };

  
    void deallocate() {
      // printf("deallocate: %x\n", device_ptr);
        if (device_ptr != nullptr) {
            cudaError_t ret = cudaFree(device_ptr);
            if (ret != cudaSuccess)
                printf("CudaDeviceMemory: device free error\n");
            device_ptr = nullptr;
            size       = 0;
        }
    };

    bool allocate(size_t size_in) {
        if (device_ptr != nullptr) {
            deallocate();
        }
        cudaError_t ret = cudaMalloc((void**)&device_ptr, size_in * sizeof(T));

        if (ret == cudaSuccess) {
            size = size_in;
        }
        else {
            size       = 0;
            device_ptr = nullptr;
        }

        return ret == cudaSuccess;
    };

    T* ptr() {
        return device_ptr;
    };

    T* operator*() {
        return device_ptr;
    };

    bool allocated() {
        return device_ptr != nullptr;
    };

    size_t get_size() {
        return size;
    };


    // zero out device memory
    bool zero() {

        cudaError_t ret = cudaMemset(device_ptr, 0, sizeof(T) * size);
        return ret == cudaSuccess;
    };

    // copy data from device to local array
    bool fetch(std::unique_ptr<T[]>& host_ptr) {
        cudaError_t ret =
            cudaMemcpy(host_ptr.get(), device_ptr, size * sizeof(T), cudaMemcpyDeviceToHost);
        return ret == cudaSuccess;
    };

    // copy data from local array to device
    bool put(std::unique_ptr<T[]>& host_ptr) {
        cudaError_t ret =
            cudaMemcpy(device_ptr, host_ptr.get(), size * sizeof(T), cudaMemcpyHostToDevice);
        return ret == cudaSuccess;
    };

private:
  void register_to_cdmm() {
    cuda_device_memory_manager & cdmm = cuda_device_memory_manager::get_instance();
    cdmm.register_mem(this);
  }
  void unregister_from_cdmm() {
    cuda_device_memory_manager & cdmm = cuda_device_memory_manager::get_instance();
    cdmm.unregister_mem(this);
  }

  
    T*     device_ptr = nullptr;
    size_t size;
};

