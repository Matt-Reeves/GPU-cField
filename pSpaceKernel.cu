#include <getCUDA.h>
#include <useful.h>
#include <grid.h>
#include <deviceFunctions.cuh>


__global__ void pSpaceKernel( int n, double nx, double ny, double dx, double dy, double x0, double y0, double renorm, double2* psi )
{
  int i;

  double x, y, pot;//,nx,di,pot;
  
//  double2 p;

  for (i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
  {
    // Start loading psi[i]... does this speed things up?
//    p = psi[i];
    // Working 1D expression
    //x = fma( DBLE(i), deltax, offsetx );
    //psi[i].x = x;
    //psi[i].y = 0.0;

    // IF THERE IS A POTENTIAL, CALCULATE IT!
    // Working 2D expression
    x = fmod( __int2double_rn(i), nx );
    y = ceil((__int2double_rn(i) - x)/nx);
    x = fma( x, dx, x0 );
    y = fma( y, dy, y0 );
    //pot = 0.0;//x*x + y*y;
    pot  = 3.0*(1+tanh((sqrt(x*x+y*y)-200.0)/1.0))/2;
    //pot += 10*exp(-(x*x+y*y)/13.1802045796452/13.1802045796452);
    pot += 20*(1+tanh((25-sqrt(x*x+y*y))/1.0))/2;
// Calculate position space terms
    psi[i] = addc( mulrc(pow(renorm,3) * modsq(psi[i]), psi[i]) , mulrc(renorm*(pot), psi[i]) );
  }

}

