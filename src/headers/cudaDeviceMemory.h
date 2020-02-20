#pragma once

#include <cstddef>
#include <cstdio>
#include <memory>

// class to manage device memory, to take care of allocation and deallocation.
template<typename T> class cuda_device_memory
{
public:
    cuda_device_memory(){

    };

    cuda_device_memory(size_t size) {
        allocate(size);
    }

    ~cuda_device_memory() {
        deallocate();
    };

    void deallocate() {
        if (device_ptr != nullptr) {
            cudaError_t ret = cudaFree(device_ptr);
            if (ret != cudaSuccess)
                printf("CudaDeviceMemory: device free error\n");
            device_ptr = nullptr;
            size       = 0;
        }
    }

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
    T*     device_ptr = nullptr;
    size_t size;
};
