#include <getCUDA.h>
#include <useful.h>
#include <grid.h>
#include <deviceFunctions.cuh>


__global__ void mulrcKernel( int n, double renorm, double2* phi)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int i;

  for (i = tid; i < n; i += blockDim.x * gridDim.x)
  { 
    phi[i] = mulrc(renorm, phi[i]);
  }

}


