
import bindtest as lib
import numpy as np

# CUDA
import pycuda.driver as cuda
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
import pycuda.driver as drv
from pycuda.compiler import SourceModule

# CFFI

import ctypes

# pybind11 version

# CUDA test
a = np.ones((4, 4), dtype=np.float64)*1.5
b = np.arange(16, dtype=np.float64).reshape(4, 4)

c = np.zeros((4, 4), dtype=np.float64)

mod = SourceModule("""
  __global__ void summy(double *a, double *b, double *c, int s_x, int s_y)
  {
    int idx = threadIdx.x + threadIdx.y*s_x;
    c[idx] = a[idx] + b[idx];
  }
  """)


#func = mod.get_function("summy")
# func(drv.In(a), drv.In(b), drv.Out(c), np.int32(
#    4), np.int32(4), block=(4, 4, 1), grid=(1, 1))
#
# print(c)


print("printing test")
lib.print_hello()

aplusb = lib.add_test(2, 4)
print(f"a + b = {aplusb}")

print("CUDA test")

a_gpu = gpuarray.to_gpu(a)
b_gpu = gpuarray.to_gpu(b)
c_gpu = gpuarray.to_gpu(c)

# lib.call_summy(ffi.cast("double *", a_gpu.ptr),
#               ffi.cast("double *", b_gpu.ptr),
#               ffi.cast("double *", c_gpu.ptr),
#               4, 4)

print(f"a_gpu pointer: {a_gpu.ptr:x}")
print(f"b_gpu pointer: {b_gpu.ptr:x}")
print(f"c_gpu pointer: {c_gpu.ptr:x}")

lib.call_summy(a_gpu.ptr,
               b_gpu.ptr,
               c_gpu.ptr,
               4, 4)


c_res = c_gpu.get()
print(c_res)
