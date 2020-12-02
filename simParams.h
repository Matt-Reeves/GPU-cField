/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_SIM_PARAMS_H
#define HAS_SIM_PARAMS_H

#include<cufft.h>
typedef struct simParams {
  double delta_t;
  double time;
  double end_time;
  double gamma;
  double mu;
  double tol;
  double renorm_p_space;
  double renorm_k_space;
  double oldmax;
  double ksqmax;
  cufftHandle plan;
} simParams;


#endif
