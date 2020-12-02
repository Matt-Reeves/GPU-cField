/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_P_SPACE_KERNEL_H
#define HAS_P_SPACE_KERNEL_H

//#include <grid.h>
//__global__ void pSpaceKernel( grid gd, double renorm, double2* psi );
//
#include<getCUDA.h>
__global__ void pSpaceKernel( int n, double nx, double ny, double dx, double dy, double x0, double y0, double renorm, double2* psi );

#endif

