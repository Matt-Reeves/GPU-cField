#include <getCUDA.h>
#include <useful.h>
#include <grid.h>
#include <deviceFunctions.cuh>


__global__ void initialProjectorKernel(int n, double nx, double ny, double dkx, double dky, double renorm, double ksqmax, double2* phi)
{
  int i;
  //int n = gd.n.x * gd.n.y * gd.n.z;
  double kx, ky;
  double pj;
  //double ky,di,ny,indx,pj;
  //nx = __int2double_rn(gd.n.x);
  //ny = DBLE(gd.n.y);
  for (i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
  { 

    // Get kx and ky
//    di = DBLE(i);
//    indx = fmod( di, nx );
    kx = fmod(__int2double_rn(i),nx);
    ky = ceil((__int2double_rn(i) - kx)/nx);

    kx = dkx * (kx - (nx * 0.5 *(1.0+tanh(BIG*(kx-0.5*(nx+1.0))))));
    ky = dky * (ky - (ny * 0.5 *(1.0+tanh(BIG*(ky-0.5*(ny+1.0))))));

    //ksq = kx*kx + ky*ky;

    // Projector
    pj = 0.5 * (1.0 - tanh(BIG*(pow(kx+LITTLE,2)+pow(ky+LITTLE,2)-ksqmax))); // NOTE: ksqmax should be passed in a little too big...
    //pj = 1.0;
    phi[i] = mulrc(renorm*pj, phi[i]);

    //phi[i].x = pj;
    //phi[i].y = 0.0;

 }

}

