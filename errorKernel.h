#ifndef HAS_ERROR_KERNEL_H
#define HAS_ERROR_KERNEL_H

#include<getCUDA.h>
__global__ void errorKernel( int n, double2** arrs, double* err, double* rho, double oldmax );

#endif

