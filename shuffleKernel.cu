#include <getCUDA.h>
#include <useful.h>
#include <deviceFunctions.cuh>
#include <cashCarp45.h>
#include <stdio.h>


template <int stage> __global__ void shuffleKernel( int n, double2** arrs)
{
  int i;
  double2 res;

  for (i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
  { 
    // Update forwards
    res = arrs[1][i];
    if (stage == 1){
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
        //arrs[7][i] = addc( arrs[7][i], mulrc(D2,res) );
        //arrs[8][i] = addc( arrs[8][i], mulrc(C2,res) );
    }
    if (stage == 3){
        arrs[3][i] = res;
        arrs[4][i] = addc( addc( arrs[4][i], mulrc(B42,arrs[2][i]) ), mulrc(B43,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D3,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C3,res) );
    }
    if (stage == 4){
        arrs[4][i] = res;
        arrs[5][i] = addc( addc( addc( arrs[5][i], mulrc(B52,arrs[2][i]) ), mulrc(B53,arrs[3][i]) ), mulrc(B54,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D4,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C4,res) );
    }
    if (stage == 5){
        arrs[6][i] = addc( addc( addc( addc( arrs[6][i], mulrc(B62,arrs[2][i]) ), mulrc(B63,arrs[3][i]) ), mulrc(B64,arrs[4][i]) ), mulrc(B65,res) );
        arrs[7][i] = addc( arrs[7][i], mulrc(D5,res) );
        //arrs[8][i] = addc( arrs[8][i], mulrc(C5,res) );
    }
    if (stage == 6){
        arrs[7][i] = addc( arrs[7][i], mulrc(D6,res) );
        arrs[8][i] = addc( arrs[8][i], mulrc(C6,res) );
    }

  }

}



void call_shuffleKernel(int stage, int num_blocks, int threads_per_block, int n, double2** arrs)
{
  if (stage == 1) {shuffleKernel<1><<< num_blocks, threads_per_block >>>( n, arrs);}
  if (stage == 2) {shuffleKernel<2><<< num_blocks, threads_per_block >>>( n, arrs);}
  if (stage == 3) {shuffleKernel<3><<< num_blocks, threads_per_block >>>( n, arrs);}
  if (stage == 4) {shuffleKernel<4><<< num_blocks, threads_per_block >>>( n, arrs);}
  if (stage == 5) {shuffleKernel<5><<< num_blocks, threads_per_block >>>( n, arrs);}
  if (stage == 6) {shuffleKernel<6><<< num_blocks, threads_per_block >>>( n, arrs);}
}

