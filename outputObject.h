/********************************************************************\
Copyright 2015 Thomas P. Billam, Durham University

This file is part of the cFieldGPU package.
\********************************************************************/

#ifndef HAS_OUTPUT_OBJECT_H
#define HAS_OUTPUT_OBJECT_H

class outputObject {
public:
  // Define an interface for outputting a simulation timestep
  // Hopefully this is sufficiently generic for many purposes
  virtual int outputByIndex( int index, grid gd, double time, double* psi ) = 0;

};

#endif

