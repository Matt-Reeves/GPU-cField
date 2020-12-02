#include <getCUDA.h>
#include <useful.h>
#include <deviceFunctions.cuh>
#include <stdio.h>


__global__ void errorKernel( int n, double2** arrs, double* err, double* rho, double oldmax)
{ 
  int i;
  //double pj;

  for (i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
  {
    
    //pj = 1.0;//rho[i];
    //err[i] = (absdiffc( arrs[8][i], arrs[7][i] )/pj) * 0.5 * (1.0 + tanh(BIG*(pj-ERR_CUT_OFF*oldmax)));
    err[i] = ( absdiffc( arrs[8][i], arrs[7][i] ) );
  }

}

