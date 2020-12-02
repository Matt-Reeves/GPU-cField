#include <getCUDA.h>
#include <useful.h>
#include <grid.h>
#include <deviceFunctions.cuh>
#include <cashCarp45.h>
#include <stdio.h>


template <int stage> __global__ void kSpaceKernel( int n, double nx, double ny, double dkx, double dky, double renorm, double2 timeconst, double ksqmax, double mu,
                                                   double2** arrs)
{
  int i;
  double kx, ky;
  double pj;
  double2 res2;//, res;

  for (i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
  { 

    kx = fmod(__int2double_rn(i),nx);
    ky = ceil((__int2double_rn(i) - kx)/nx);

    kx = dkx * (kx - (nx * 0.5 *(1.0+tanh(BIG*(kx-0.5*(nx+1.0))))));
    ky = dky * (ky - (ny * 0.5 *(1.0+tanh(BIG*(ky-0.5*(ny+1.0))))));

    pj = (0.5 * (1.0 - tanh(BIG*((pow(kx+LITTLE,2)+pow(ky+LITTLE,2))-ksqmax)))); // NOTE: ksqmax should be passed in a little too big...
    res2 = mulrc( pj, timeconst);
    //res2 = timeconst;

    // Input variable "timeconst" is expected to have gamma already in it...
//    if (stage == 1)
//    {
//      arrs[1][i] = mulc( res2 , addc( mulrc(renorm,arrs[1][i]) , mulrc(0.5*(pow(kx,2)+pow(ky,2))-mu,arrs[0][i]) ) );
//    }
//    else
//    {
      arrs[1][i] = mulc( res2 , addc( mulrc(renorm,arrs[1][i]) , mulrc(0.5*(pow(kx,2)+pow(ky,2))-mu,arrs[stage][i]) ) );
//    }
/*
    // Update forwards
    if (stage == 1){
        //arrs[1][i] = res; // weirdly, this step turns out not to be necessary!
        arrs[2][i] = addc( arrs[0][i], mulrc(B21,res) );
        arrs[3][i] = addc( arrs[0][i], mulrc(B31,res) );
        arrs[4][i] = addc( arrs[0][i], mulrc(B41,res) );
        arrs[5][i] = addc( arrs[0][i], mulrc(B51,res) );
        arrs[6][i] = addc( arrs[0][i], mulrc(B61,res) );
        arrs[7][i] = addc( arrs[0][i], mulrc(D1,res) );
        arrs[8][i] = addc( arrs[0][i], mulrc(C1,res) );
    }
    if (stage == 2){
        arrs[2][i] = res;
        arrs[3][i] = addc( arrs[3][i], mulrc(B32,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D2,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C2,res) );
    }
    if (stage == 3){
        arrs[3][i] = res;
        arrs[4][i] = addc( addc( arrs[3][i], mulrc(B42,arrs[2][i]) ), mulrc(B43,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D3,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C3,res) );
    }
    if (stage == 4){
        arrs[4][i] = res;
        arrs[5][i] = addc( addc( addc( arrs[4][i], mulrc(B52,arrs[2][i]) ), mulrc(B53,arrs[3][i]) ), mulrc(B54,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D4,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C4,res) );
    }
    if (stage == 5){
        arrs[5][i] = res;
        arrs[6][i] = addc( addc( addc( addc( arrs[5][i], mulrc(B62,arrs[2][i]) ), mulrc(B63,arrs[3][i]) ), mulrc(B64,arrs[4][i]) ), mulrc(B65,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D5,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C5,res) );
    }
    if (stage == 6){
        // Get the 5th-order accurate solution
        arrs[8][i] = addc( arrs[8][i], mulrc(C6,res) );
   
        //res2 = addc( arrs[8][i], mulrc(C6,res) );
        //arrs[8][i] = res2;
        //pj = sqrt(modsq(res2));
        //err[i] = (absdiffc( res2, addc( arrs[7][i], mulrc(D6,res) ) )/pj) * 0.5 * (1.0 + tanh(BIG*(pj-ERR_CUT_OFF*oldmax))) ;
        //rho[i] = pj;
   
        // Get the 4th-order solution (not usually necessary)
        arrs[7][i] = addc( arrs[7][i], mulrc(D6,res) );
    }
*/

  }
}

void call_kSpaceKernel(int stage, int num_blocks, int threads_per_block, int n, double nx, double ny, double dkx, double dky, double renorm, double2 timeconst, double ksqmax, double mu, double2** arrs)
{
  if (stage == 1) {kSpaceKernel<0><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
  if (stage == 2) {kSpaceKernel<2><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
  if (stage == 3) {kSpaceKernel<3><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
  if (stage == 4) {kSpaceKernel<4><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
  if (stage == 5) {kSpaceKernel<5><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
  if (stage == 6) {kSpaceKernel<6><<< num_blocks, threads_per_block >>>( n, nx, ny, dkx, dky, renorm, timeconst, ksqmax, mu, arrs);}
}



