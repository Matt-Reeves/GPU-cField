/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_SHUFFLE_KERNEL_H
#define HAS_SHUFFLE_KERNEL_H

#include <getCUDA.h>

template <int stage> __global__ void shuffleKernel( int n, double2** arrs);

void call_shuffleKernel(int stage, int num_blocks, int threads_per_block, int n, double2** arrs);

#endif

