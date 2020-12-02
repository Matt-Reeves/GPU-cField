#include <getCUDA.h>

__global__ void copyKernel( int n, double2* a, double2* b )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int i;
  for (i = tid; i < n; i += blockDim.x * gridDim.x)
  {
    b[i] = a[i];
  }
}

