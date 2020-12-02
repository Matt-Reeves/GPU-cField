/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_PGP_EVOLVER
#define HAS_PGP_EVOLVER

#include <math.h>
#include <getCUDA.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <stdio.h>
#include <stdlib.h>
#include <getCUDA.h>
#include <pSpaceKernel.h>
#include <kSpaceKernel.h>
#include <initialProjectorKernel.h>
#include <mulrcKernel.h>
#include <grid.h>
#include <useful.h>
#include <copyKernel.h>
#include <absArrayKernel.h>
#include <errorKernel.h>
#include <shuffleKernel.h>
#include <cufft.h>
#include <outputObject.h>
#include <time.h>

class PGPEvolver
{
  public:
  // Constructor
  PGPEvolver( grid gd_in, double2 *psi, double time_in, outputObject* obj, double gamma_in, double delta_t_in, int zero_index_in, int threads_per_block_in, int num_blocks_in );

  void evolveSteps( int nsteps, double output_delta_t, double max_seconds );

  void evolveToTime( double end_time, int& gc, int& bc );

  void doOutput(int i);

  ~PGPEvolver( void );


  private:
  static const double mu = 1.0;
  static const double tol = 1.0E-6;
  static const double safety = 0.9;
  

 // bool ready_to_evolve; //??????

  double2 *psi;
  outputObject *op;

  int zero_index;
  grid gd;
  double t; 
  int n;
  double renorm_pspace;
  double renorm_kspace;
  double ksqmax; // ** NOTE: THIS ONLY WORKS FOR SQUARE GRID
  double oldmax;
  double gamma;
  double2 *a, **arrs, **arrs_host;
  double *err, *rho;
  thrust::device_ptr<double> err_thrust;
  thrust::device_ptr<double> rho_thrust;
  cufftHandle plan;
  double delta_t;  
  double fnx, fny, fdx, fdy, fx0, fy0, fdkx, fdky;
  int threads_per_block;
  int num_blocks; 

};

#endif

