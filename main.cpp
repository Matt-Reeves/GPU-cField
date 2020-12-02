#include <PGPEvolver.h>
#include <HDF5io.h>
#include <stdio.h>
#include <getCard.h>
#include <describeVersion.h>

int main(int argc, char* argv[])
{
  // Print version string
  printf("This is cFieldGPU, Version %s\n", VERSION);

  // Which card?
  if ( argc < (1 + 5) )
  {
    printf("Not enough arguments (need cardNum, gamma, delta_t, num_steps, max_hours)\n");
    return -1;
  }

  // Parse command line args
  int cardNum;
  sscanf(argv[1],"%d",&cardNum);
  printf("Read in device id %d\n",cardNum);

  double gamma;
  sscanf(argv[2],"%lf",&gamma);
  printf("Read in gamma %f\n",gamma);

  double delta_t;
  sscanf(argv[3],"%lf",&delta_t);
  printf("Read in delta_t %f\n",delta_t);

  int num_steps;
  sscanf(argv[4],"%d",&num_steps);
  printf("Read in num_steps %d\n",num_steps);

  double max_hours;
  sscanf(argv[5],"%lf",&max_hours);
  printf("Read in max_hours %f\n",max_hours);

  double *phi;
  grid gd;
  char filename[255];
  int dim, err;
  HDF5io h;

  if (!sprintf(filename,"0.h5"))
  {
     printf("Error writing filename to string\n");
     return -1;
  }

  dim = h.readAttributes(filename,&gd);
  if (dim < 0)
  {
    printf("No good dimensions in file!\n");
    return -1;
  }
  printf("Got %d-dimensional grid\n", dim);
  printf("gd.n.x: %d\n",gd.n.x);
  printf("gd.delta.x: %g\n",gd.delta.x);
  printf("gd.offset.x: %g\n",gd.offset.x);
  printf("gd.n.y: %d\n",gd.n.y);
  printf("gd.delta.y: %g\n",gd.delta.y);
  printf("gd.offset.y: %g\n",gd.offset.y);
  printf("gd.n.z: %d\n",gd.n.z);
  printf("gd.delta.z: %g\n",gd.delta.z);
  printf("gd.offset.z: %g\n",gd.offset.z);


printf("Hello World!!!!\n");
  // Set up rest of grid
  gd.dk.x = 2.0*PI/(gd.delta.x *((double)gd.n.x));
  gd.dk.y = 2.0*PI/(gd.delta.y *((double)gd.n.y));

  // Allocate array
  phi = (double *) malloc(2*gd.n.x*gd.n.y*gd.n.z*sizeof(double));

  err = h.readHDF5Wavefunction2D(filename, phi, gd);
  if (err < 0)
  {
    printf("Read failed!\n");
    return -1;
  }

  // Get a GPU to run on
  getCard(cardNum);
  // Initialize an evolver
  PGPEvolver p = PGPEvolver( gd, (double2*) phi, 0.0, &h, gamma, 0.00001, 1, 1024, 1024 ); // Note cast from double to double2

  p.evolveSteps( num_steps, delta_t, max_hours*3600.0);

  // Output
/*  if (!sprintf(filename,"end.h5"))
  {
     printf("Error writing filename to string\n");
     return -1;
  }
  err = h.writeHDF5Wavefunction2D(filename, phi, gd);
  if (err < 0)
  {
    printf("Write failed!\n");
    return -1;
  }
*/

  free( phi );
  return 0;
}

