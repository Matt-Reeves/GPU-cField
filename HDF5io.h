#ifndef HAS_HDF5_IO_H
#define HAS_HDF5_IO_H
#include <grid.h>
#include <outputObject.h>
#include "hdf5.h"

class HDF5io: public outputObject 
{
  public:
  // Overrride interface specified by outputObject
  int outputByIndex( int index, grid gd, double time, double* psi );

  // Provide some other stuff
  int outputByName( char* name, grid gd, double time, double* psi );
  int inputByName( char* name, grid gd, double time, double* psi );
  int inputByIndex( int index, grid gd, double time, double* psi );

  int readAttributes(char* filename, grid* gd);
  int readHDF5Wavefunction2D(char* filename, double* phi, grid gd);
  int writeHDF5Wavefunction2D(char* filename, double* phi, grid gd);

};

#endif
