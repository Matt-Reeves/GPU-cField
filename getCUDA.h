#ifndef HAS_CUDA
#define HAS_CUDA

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#define cudaErr(result) {cudaErrorCheck(result, __FILE__, __LINE__);}

inline void cudaErrorCheck(cudaError_t result, char *file, int line, bool abort=true)
{
   if (result != cudaSuccess) 
   {
      fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(result), file, line);
      if (abort) 
      {
        exit(result);
      }
   }
}

#endif

