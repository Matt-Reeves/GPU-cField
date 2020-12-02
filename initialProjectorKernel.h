#ifndef HAS_INITIAL_PROJECTOR_KERNEL_H
#define HAS_INITIAL_PROJECTOR_KERNEL_H

__global__ void initialProjectorKernel(int n, double nx, double ny, double dkx, double dky, double renorm, double ksqmax, double2* phi);

#endif

