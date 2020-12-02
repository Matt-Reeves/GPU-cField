#include <getCUDA.h>
#include <iostream>
#include <stdio.h>
#include <getCard.h>

void getCard(int cardNum)
{
  int numCards, gotCardNum;
  size_t freeMem, totalMem, usedMem;
  cudaDeviceProp cardInfo;
  
  // Get card
  cudaErr( cudaGetDeviceCount(&numCards) );
  printf("Found %d GPU cards\n",numCards);
  printf("Attempting to select card id %d\n", cardNum);
  cudaErr( cudaSetDevice(cardNum) );
  cudaErr( cudaFree(0) );
  cudaErr( cudaGetDevice(&gotCardNum) );
  printf("Selected card id %d\n", gotCardNum);
  cudaErr( cudaGetDeviceProperties(&cardInfo, gotCardNum) );
  printf("Card type: %s\n", cardInfo.name);
  printf("Total mem: %lu MB\n", cardInfo.totalGlobalMem/1024/1024);
  printf("SMX units: %d\n", cardInfo.multiProcessorCount);
  printf("Max block: %d threads\n",cardInfo.maxThreadsPerBlock); 
  printf("Share mem: %lu KB per block\n", cardInfo.sharedMemPerBlock/1024);
  printf("Registers: %d per block\n", cardInfo.regsPerBlock);
  printf("Warp size: %d\n", cardInfo.warpSize);
  printf("Clock MHz: %d\n", cardInfo.clockRate/1000);
  // Summarize memory information
  cudaErr( cudaMemGetInfo( &freeMem, &totalMem ) ) ;
  usedMem = totalMem - freeMem;
  printf(" Free mem: %lu MB\n", freeMem/1024/1024);
  printf(" Used mem: %lu MB\n", usedMem/1024/1024);

  return;
}
 

