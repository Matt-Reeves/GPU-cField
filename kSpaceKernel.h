#ifndef HAS_K_SPACE_KERNEL_H
#define HAS_K_SPACE_KERNEL_H

#include <getCUDA.h>
template <int stage> __global__ void kSpaceKernel( int n, double nx, double ny, double dkx, double dky, double renorm, double2 timeconst, double ksqmax, double mu, double2** arrs);

void call_kSpaceKernel(int stage, int num_blocks, int threads_per_block, int n, double nx, double ny, double dkx, double dky, double renorm, double2 timeconst, double ksqmax, double mu, double2** arrs);

#endif

