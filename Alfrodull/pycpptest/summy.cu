
#include <cstdio>

#include <pybind11/pybind11.h>


void print_hello() {
  printf("Hello World!\n");
}

int add_test(int a, int b)
{
  return a + b;
}



__global__ void summy(double *a, double *b, double *c, int s_x, int s_y) {
  int idx = threadIdx.x + threadIdx.y * s_x;
  c[idx] = a[idx] + b[idx];
}

void call_summy(long a_, long b_, long c_, int s_x, int s_y) {
  double * a = (double *)a_;
  double * b = (double *)b_;
  double * c = (double *)c_;
  dim3 threadsPerBlock(1, 1);
  dim3 numBlocks(s_x, s_y, 1);
  printf("a: %p, b: %p, c: %p\n", a, b, c);
  summy<<<threadsPerBlock, numBlocks>>>(a, b, c, s_x, s_y);

  cudaDeviceSynchronize();
}

PYBIND11_MODULE(bindtest, m) {
    m.doc() = "pybind11 test for cuda wrapping plugin"; // optional module docstring

    m.def("add_test", &add_test, "A function which adds two numbers");
    m.def("print_hello", &print_hello, "A function that prints hello");
    m.def("call_summy", &call_summy, "A function that calls a cuda sum function");
}
