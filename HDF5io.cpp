#include <stdio.h>
#include <HDF5io.h>
#include <stdlib.h>
#include "hdf5.h"
#include <grid.h>

void writeAttrInt(hid_t gid, const char* name, int val)
{
  hid_t aid, aspace;
  herr_t err;

  aspace = H5Screate(H5S_SCALAR);
  aid = H5Acreate(gid,name,H5T_NATIVE_INT,aspace, H5P_DEFAULT, H5P_DEFAULT);
  err = H5Awrite(aid, H5T_NATIVE_INT, &val);
  err = H5Aclose(aid);
  err = H5Sclose(aspace);
}

void writeAttrDouble(hid_t gid, const char* name, double val)
{
  hid_t aid, aspace;
  herr_t err;

  aspace = H5Screate(H5S_SCALAR);
  aid = H5Acreate(gid,name,H5T_NATIVE_DOUBLE,aspace, H5P_DEFAULT, H5P_DEFAULT);
  err = H5Awrite(aid, H5T_NATIVE_DOUBLE, &val);
  err = H5Aclose(aid);
  err = H5Sclose(aspace);
}


void stopReading(hid_t gid, hid_t fid)
{
  herr_t err;
  // Closing down
  err = H5Gclose(gid);
  if (err < 0)
  {
    printf("Error closing group 1\n");
  }
  err = H5Fclose(fid);
  if (err < 0)
  {
    printf("Error closing h5 file\n");
  }
}

int HDF5io::readAttributes(char* filename, grid* gd)
{
  hid_t fid, gid, aid;
  herr_t err;

  printf("%s\n",filename);
  (*gd).n.x = 1;
  (*gd).n.y = 1;
  (*gd).n.z = 1;


  fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0)
  {
    printf("Error opening h5 file %s\n",filename);
    return -1;
  }

  gid = H5Gopen( fid, "/1", H5P_DEFAULT );
  if (gid < 0)
  {
    printf("Error opening group 1");
    return -1;
  } 

  // Try to get all the x attributes
  aid = H5Aopen( gid, "nx", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute nx\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Aread( aid, H5T_NATIVE_INT, &(*gd).n.x );
  if (err < 0)
  {
    printf("Error reading attribute nx\n");
    stopReading(gid,fid);
    return -1;
  }

  aid = H5Aopen( gid, "dx", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute dx\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).delta.x );
  if (err < 0)
  {
    printf("Error reading attribute dx\n");
    stopReading(gid,fid);
    return -1;
  }

  aid = H5Aopen( gid, "x0", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute x0\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).offset.x );
  if (err < 0)
  {
    printf("Error reading attribute x0\n");
    stopReading(gid,fid);
    return -1;
  }


  // Try to get all the y attributes
  aid = H5Aopen( gid, "ny", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute ny\n");
    stopReading(gid,fid);
    return 1;
  }

  err = H5Aread( aid, H5T_NATIVE_INT, &(*gd).n.y );
  if (err < 0)
  {
    printf("Error reading attribute ny\n");
    stopReading(gid,fid);
    return 1;
  }

  aid = H5Aopen( gid, "dy", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute dy\n");
    stopReading(gid,fid);
    return 1;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).delta.y );
  if (err < 0)
  {
    printf("Error reading attribute dy\n");
    stopReading(gid,fid);
    return 1;
  }

  aid = H5Aopen( gid, "y0", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute y0\n");
    stopReading(gid,fid);
    return 1;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).offset.y );
  if (err < 0)
  {
    printf("Error reading attribute y0\n");
    stopReading(gid,fid);
    return 1;
  }



  // Try to get all the z attributes
  aid = H5Aopen( gid, "nz", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute nz\n");
    stopReading(gid,fid);
    return 2;
  }

  err = H5Aread( aid, H5T_NATIVE_INT, &(*gd).n.z );
  if (err < 0)
  {
    printf("Error reading attribute nz\n");
    stopReading(gid,fid);
    return 2;
  }

  aid = H5Aopen( gid, "dz", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute dz\n");
    stopReading(gid,fid);
    return 2;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).delta.z );
  if (err < 0)
  {
    printf("Error reading attribute dz\n");
    stopReading(gid,fid);
    return 2;
  }

  aid = H5Aopen( gid, "z0", H5P_DEFAULT );
  if (aid < 0)
  {
    printf("Error opening attribute z0\n");
    stopReading(gid,fid);
    return 2;
  }

  err = H5Aread( aid, H5T_NATIVE_DOUBLE, &(*gd).offset.z );
  if (err < 0)
  {
    printf("Error reading attribute z0\n");
    stopReading(gid,fid);
    return 2;
  }


  stopReading(gid,fid);
  return 3;
  
}


int HDF5io::readHDF5Wavefunction2D(char* filename, double* phi, grid gd)
{
  hid_t fid, gid, memspace_id, dataset_id;
  herr_t err;
  hsize_t dims[2], stride[1], off_imag[1], off_real[1], mem_dims[1], length[1], block[1];

  dims[0] = gd.n.x;
  dims[1] = gd.n.y;
  mem_dims[0] = 2*gd.n.x*gd.n.y;
  stride[0] = 2;
  length[0] = gd.n.x*gd.n.y;
  off_real[0] = 0;
  off_imag[0] = 1;
  block[0] = 1;

  printf("%s\n",filename);

  /* Open file and group */
  fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0)
  {
    printf("Error opening h5 file %s\n",filename);
    return -1;
  }
  gid = H5Gopen( fid, "/1", H5P_DEFAULT );
  if (gid < 0)
  {
    printf("Error opening group 1\n");
    return -1;
  } 
  /* ------------------ */

  memspace_id = H5Screate_simple( 1, mem_dims, mem_dims );
  if (memspace_id < 0)
  {
    printf("Failed to create memory dataspace...\n");
    stopReading(gid,fid);
  }

  /* IMAGINARY DATA */
  dataset_id = H5Dopen(gid, "phiI", H5P_DEFAULT);
  if (dataset_id < 0)
  {
    printf("Error opening imaginary dataset in file\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, off_imag, stride, length, block);
  if (err < 0)
  {
    printf("Error selecting imaginary hyperslab in memory\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, H5S_ALL, H5P_DEFAULT , phi);
  if (err < 0)
  {
    printf("Error reading imaginary hyperslab\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dclose(dataset_id);
  if (err < 0)
  {
    printf("Couldn't close imaginary dataset in file...\n");
  }
  /* ---------------------- */

  /* REAL DATA */
  dataset_id = H5Dopen(gid, "phiR", H5P_DEFAULT);
  if (dataset_id < 0)
  {
    printf("Error opening real dataset in file\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, off_real, stride, length, block);
  if (err < 0)
  {
    printf("Error selecting real hyperslab in memory\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, H5S_ALL, H5P_DEFAULT , phi);
  if (err < 0)
  {
    printf("Error reading real hyperslab\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dclose(dataset_id);
  if (err < 0)
  {
    printf("Error closing real dataset\n");
  }
  /* ---------------------- */



  err = H5Sclose(memspace_id);
  if (err < 0)
  {
    printf("Couldn't close memory dataspace...\n");
  }

  /* Close file and group */
  stopReading(gid,fid);
  /* ------------------- */
  return 0;
}


int HDF5io::writeHDF5Wavefunction2D(char* filename, double* phi, grid gd)
{
  hid_t fid, gid, filespace_id, dataset_id, memspace_id;
  herr_t err;
  hsize_t dims[2], stride[1], off_imag[1], off_real[1], mem_dims[1], length[1], block[1];

  dims[0] = gd.n.x;
  dims[1] = gd.n.y;
  mem_dims[0] = 2*gd.n.x*gd.n.y;
  stride[0] = 2;
  length[0] = gd.n.x*gd.n.y;
  off_real[0] = 0;
  off_imag[0] = 1;
  block[0] = 1;

  printf("%s\n",filename);

  /* Open file and group */
  fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0)
  {
    printf("Error creating h5 file %s\n",filename);
    return -1;
  }
  gid = H5Gcreate( fid, "/1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  if (gid < 0)
  {
    printf("Error creating group 1\n");
    return -1;
  } 
  /* ------------------ */

  memspace_id = H5Screate_simple( 1, mem_dims, mem_dims );
  if (memspace_id < 0)
  {
    printf("Failed to create memory dataspace...\n");
    stopReading(gid,fid);
  }

  filespace_id = H5Screate_simple( 2, dims, dims );
  if (filespace_id < 0)
  {
    printf("Failed to create file dataspace...\n");
    stopReading(gid,fid);
  }

  /* IMAGINARY DATA */
  dataset_id = H5Dcreate(gid, "phiI", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0)
  {
    printf("Error creating imaginary dataset in file\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, off_imag, stride, length, block);
  if (err < 0)
  {
    printf("Error selecting imaginary hyperslab in memory\n");
    stopReading(gid,fid);
    return -1;
  }
  err = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT , phi);
  if (err < 0)
  {
    printf("Error writing imaginary hyperslab\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dclose(dataset_id);
  if (err < 0)
  {
    printf("Couldn't close imaginary dataset in file...\n");
  }
  /* ---------------------- */

  /* REAL DATA */
  dataset_id = H5Dcreate(gid, "phiR", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0)
  {
    printf("Error creating real dataset in file\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, off_real, stride, length, block);
  if (err < 0)
  {
    printf("Error selecting real hyperslab in memory\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, H5P_DEFAULT , phi);
  if (err < 0)
  {
    printf("Error writing real hyperslab\n");
    stopReading(gid,fid);
    return -1;
  }

  err = H5Dclose(dataset_id);
  if (err < 0)
  {
    printf("Error closing real dataset\n");
  }
  /* ---------------------- */



  err = H5Sclose(memspace_id);
  if (err < 0)
  {
    printf("Couldn't close memory dataspace...\n");
  }
  err = H5Sclose(filespace_id);
  if (err < 0)
  {
    printf("Couldn't close filespace...\n");
  }

  
  /* WRITE ATTRIBUTES */
  writeAttrInt(gid, "nx", gd.n.x);
  writeAttrInt(gid, "ny", gd.n.y);
  writeAttrDouble(gid, "dx", gd.delta.x);
  writeAttrDouble(gid, "dy", gd.delta.y);
  writeAttrDouble(gid, "x0", gd.offset.x);
  writeAttrDouble(gid, "y0", gd.offset.y);
  writeAttrDouble(gid, "t", 0.0);


  /* Close file and group */
  stopReading(gid,fid);
  /* ------------------- */
  return 0;
}


int HDF5io::outputByName( char* name, grid gd, double time, double* psi ) {
//std::cout << "Called outputByName: " << name << std::endl;
  return this->writeHDF5Wavefunction2D( name, (double*) psi, gd );

}

int HDF5io::inputByName( char* name, grid gd, double time, double* psi ) {
//std::cout << "Called inputByName: " << name << std::endl;
  return this->readHDF5Wavefunction2D( name, (double*) psi, gd);

}

int HDF5io::outputByIndex( int index, grid gd, double time, double* psi )
{
  char str[255];
  if (!sprintf(str,"%d.h5",index))
  {
    printf("Error writing filename to string\n");
    return -1;
  }

  return this->outputByName( str, gd, time, psi );
}

int HDF5io::inputByIndex( int index, grid gd, double time, double* psi ) 
{
  char str[255];
  if (!sprintf(str,"%d.h5",index)) 
  {
    printf("Error writing filename to string\n");
    return -1;
  }

  return this->inputByName( str, gd, time, psi );
}


