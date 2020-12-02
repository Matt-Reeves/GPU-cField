/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_MULRC_KERNEL_H
#define HAS_MULRC_KERNEL_H

#include <grid.h>
__global__ void mulrcKernel( int n, double renorm, double2* psi );

#endif


