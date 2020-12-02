#include <getCUDA.h>
//#include <deviceFunctions.cuh>

__global__ void absArrayKernel( int n, double2* a, double* b )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int i;
  for (i = tid; i < n; i += blockDim.x * gridDim.x)
  {
    b[i] = hypot(a[i].x,a[i].y);//sqrt(modsq(a[i]));
  }
}

