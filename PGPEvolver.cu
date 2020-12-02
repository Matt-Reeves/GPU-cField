#include <PGPEvolver.h>

PGPEvolver::PGPEvolver( grid gd_in, double2 *psi_in, double time_in, outputObject* obj, double gamma_in, double delta_t_in, int zero_index_in, int threads_per_block_in, int num_blocks_in )
{
  // FIRST, CHECK THAT WE HAVE ONLY TWO DIMENSIONS SO THAT THE PRESENT CODE CAN HANDLE IT
  int i;

  psi = psi_in;
  op = obj;

  // Set useful stuff
  gd = gd_in;
  n = gd.n.x * gd.n.y * gd.n.z;
  printf("%d\n",n);
  renorm_pspace = gd.dk.x*gd.dk.y / (2.0*PI);
  renorm_kspace = gd.delta.x*gd.delta.y / (2.0*PI);
  ksqmax = (0.5*PI/gd.delta.x)*(0.5*PI/gd.delta.x); // ** NOTE: THIS ONLY WORKS FOR SQUARE GRID
  gamma = gamma_in;
  delta_t = delta_t_in;
  t = time_in;
  zero_index = zero_index_in;
  threads_per_block = threads_per_block_in;
  num_blocks = num_blocks_in;

  // Doubles for grid computations
  fnx = (double) gd.n.x;
  fny = (double) gd.n.y;
  fdx = (double) gd.delta.x;
  fdy = (double) gd.delta.y;
  fx0 = (double) gd.offset.x;
  fy0 = (double) gd.offset.y;
  fdkx = (double) gd.dk.x;
  fdky = (double) gd.dk.y;

  // Try to allocate all the necessary memory space
  
  // Allocate one massive complex array to hold all vectors
  cudaErr( cudaMalloc( (void**) &a, sizeof( double2 ) * n * 9 ) );

  // Set up an array of pointers to the nine individual elements
  arrs_host = (double2**) malloc(9 * sizeof(double2*));
  if (arrs_host == NULL)
  {
    printf("Unable to allocate arrs_host!\n");
  }
  for (i=0; i<9; i++)
  {
    arrs_host[i] = a + i*n;
  }

  // Copy that array of pointers to the device
  cudaErr( cudaMalloc( (void**) &arrs, sizeof( double2* ) * 9 ) );
  cudaErr( cudaMemcpy(arrs, arrs_host, 9 * sizeof(double2*),cudaMemcpyHostToDevice) );

  // Allocate double arrays for error and absval
  cudaErr( cudaMalloc( (void**) &err, sizeof( double ) * n) );
  cudaErr( cudaMalloc( (void**) &rho, sizeof( double ) * n) );

  // Allocate thrust pointers to err and rho
  err_thrust = thrust::device_pointer_cast(err);
  rho_thrust = thrust::device_pointer_cast(rho);

  // Plan ffts
  cufftResult plan_err = cufftPlan2d(&plan, gd.n.x, gd.n.y, CUFFT_Z2Z);
  if (plan_err != CUFFT_SUCCESS) 
  {
    printf("Problem planning transforms...\n");
  }

  // Put field on device!!!!!!!
  cudaErr( cudaMemcpy( arrs_host[0], psi, n*sizeof(double2),cudaMemcpyHostToDevice ) );

}


void PGPEvolver::evolveSteps( int nsteps, double output_delta_t, double max_seconds )
{
  int i, gc, bc;
  time_t total_start, step_start, now;
  double total_elapsed, step_elapsed;

  // Get starting time
  total_start = time(NULL);

  for (i=0;i<nsteps;i++)
  {
    gc = 0;
    bc = 0;
    step_start = time(NULL);
    this->evolveToTime( t + output_delta_t, gc, bc );
    printf("Just finished step %d: took %d good steps, %d bad steps\n", i, gc, bc);
    this->doOutput(i);

    // Check time remaining
    now = time(NULL);
    step_elapsed = difftime(now, step_start);
    total_elapsed = difftime(now, total_start);
    printf("Last step took %g seconds (total time elapsed: %g seconds)\n", step_elapsed, total_elapsed);
    if ((max_seconds - total_elapsed) < 1.5 * step_elapsed)
    {
       // Taking another step could be risky...
      printf("There are %g seconds remaining and last step took %g secs, exiting gracefully now...\n", max_seconds-total_elapsed, step_elapsed);
      break;
    }
  }
}

void PGPEvolver::doOutput(int i)
{
  // Move array back to the host
  cudaErr( cudaMemcpy( psi, arrs_host[0], n*sizeof(double2),cudaMemcpyDeviceToHost ) );

  // Call output object
  // Check op isn't null first?
  op->outputByIndex( i+zero_index, gd, t, (double*) psi );

}


PGPEvolver::~PGPEvolver( void )
{
  
  cudaErr( cudaFree( a ) );
  cudaErr( cudaFree( arrs ) );
  cudaErr( cudaFree( err ) ); 
  cudaErr( cudaFree( rho ) );
  free( arrs_host );

}


void PGPEvolver::evolveToTime( double end_time, int& goodcount, int& failcount )
{
  cufftResult trans_err;
  int i;
  double maxerr, saved_delta_t, scale;
  bool can_has_exit;
  double2 tconst; // ** in fortran have -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t
  tconst.x = -gamma*delta_t;
  tconst.y = -delta_t;


  // Initial fft (and renorm / projection!)
  trans_err = cufftExecZ2Z(plan, arrs_host[0], arrs_host[0], CUFFT_FORWARD);
  if (trans_err != CUFFT_SUCCESS)
  {
    printf("Error on intro transform: %d\n",trans_err);
    return;
  }
  cudaErr( cudaPeekAtLastError() );
  cudaErr( cudaDeviceSynchronize() );
  #ifdef TIMING_MUCH
  cudaEvent_t start, stop;
  float kerneltime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  #endif
  initialProjectorKernel<<< num_blocks, threads_per_block >>>( n, fnx, fny, fdkx, fdky, renorm_kspace, ksqmax, arrs_host[0] );
  #ifdef TIMING_MUCH
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kerneltime, start, stop);
  printf("Elapsed time for initial projector kernel: %g ms\n",kerneltime);
  #endif
  cudaErr( cudaPeekAtLastError() );
  cudaErr( cudaDeviceSynchronize() );

  // Get oldmax using absArrayKernel and thrust::reduce
  #ifdef TIMING_MUCH
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  #endif
  absArrayKernel<<< num_blocks, threads_per_block >>>( n, arrs_host[0], rho );
  #ifdef TIMING_MUCH
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kerneltime, start, stop);
  printf("Elapsed time for abs array kernel: %g ms\n",kerneltime);
  #endif
  cudaErr( cudaPeekAtLastError() );
  cudaErr( cudaDeviceSynchronize() );
  oldmax = thrust::reduce(rho_thrust, rho_thrust + n, 0.0,thrust::maximum<double>());
  //printf("Oldmax: %g\n",oldmax);
  
  // Evolution loop
  #ifdef TIMING_MUCH
  int loopcount = 0;
  float timearray[7][5];
  int ii,jj;
  for (ii=0; ii<5; ii++)
  {
    for (jj=0; jj<7; jj++)
    {
      timearray[jj][ii] = 0.0;
    }
  }
  #endif
  while (true)
  {
    #ifdef TIMING_MUCH
    loopcount++;
    if (loopcount > 100)
    {
      break;
    }
    #endif
    can_has_exit = false;
    if ((end_time -t) < delta_t)
    {
      saved_delta_t = delta_t;
      delta_t = end_time - t;
      tconst.x = -gamma*delta_t;
      tconst.y = -delta_t;
      can_has_exit = true;
    }
  
    // ** LOOP OVER STAGES **
    for (i=1; i<7; i++)
    {
      // First transform to position
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      if (i==1)
      {
        trans_err = cufftExecZ2Z(plan, arrs_host[0], arrs_host[1], CUFFT_INVERSE);
      }
      else
      {
        trans_err = cufftExecZ2Z(plan, arrs_host[i], arrs_host[1], CUFFT_INVERSE);
      }
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Transform to position space %d took %g ms\n",i,kerneltime);
      timearray[i-1][0] +=kerneltime;
      #endif
      if (trans_err != CUFFT_SUCCESS)
      {
        printf("Error on transform to p space stage %d: %d",i,trans_err);
        return;
      }
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );
    
      // Execute pspace kernel
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      pSpaceKernel<<< num_blocks, threads_per_block >>>( n, fnx, fny, fdx, fdy, fx0, fy0, renorm_pspace, arrs_host[1] );
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Elapsed time for pSpaceKernel %d: %g ms\n",i,kerneltime);
      timearray[i-1][1] +=kerneltime;
      #endif
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );
    
      // First transform back to momentum
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      trans_err = cufftExecZ2Z(plan, arrs_host[1], arrs_host[1], CUFFT_FORWARD);
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Transform to momentum space %d took %g ms\n",i,kerneltime);
      timearray[i-1][2] +=kerneltime;
      #endif
      if (trans_err != CUFFT_SUCCESS)
      {
        printf("Error on transform to k space stage %d: %d",i,trans_err);
        return;
      }
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );
    
      // Perform kspace kernel of this stage
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      call_kSpaceKernel(i, num_blocks, threads_per_block, n, fnx, fny, fdkx, fdky, renorm_kspace, tconst, ksqmax, mu, arrs);
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Elapsed time for kSpaceKernel %d: %g ms\n",i,kerneltime);
      timearray[i-1][3] +=kerneltime;
      #endif
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );

      // Perform shuffle kernel of this stage
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      //template <int stage> __global__ void shuffleKernel( int n, double2** arrs);
      call_shuffleKernel( i, num_blocks, threads_per_block, n, arrs );
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Elapsed time for shuffleKernel %d: %g ms\n",i,kerneltime);
      timearray[i-1][4] +=kerneltime;
      #endif
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );
    }
    // Put mod of wavefunction in rho
    #ifdef TIMING_MUCH
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    #endif
    absArrayKernel<<< num_blocks, threads_per_block >>>( n, arrs_host[8], rho );
    #ifdef TIMING_MUCH
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&kerneltime, start, stop);
    //printf("Elapsed time for abs array kernel: %g ms\n",kerneltime);
    timearray[6][0] +=kerneltime;
    #endif
    cudaErr( cudaPeekAtLastError() );
    cudaErr( cudaDeviceSynchronize() );
    // Evaluate relative error array
    //printf("Evluating relative error with oldmax = %g\n", oldmax);
    #ifdef TIMING_MUCH
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    #endif
    errorKernel<<< num_blocks, threads_per_block >>>( n, arrs, err, rho, oldmax );
    #ifdef TIMING_MUCH
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&kerneltime, start, stop);
    //printf("Elapsed time for error kernel: %g ms\n",kerneltime);
    timearray[6][1] +=kerneltime;
    #endif
    cudaErr( cudaPeekAtLastError() );
    cudaErr( cudaDeviceSynchronize() );
    // Evaluate maximum relative error
    maxerr = thrust::reduce(err_thrust, err_thrust + n, 0.0,thrust::maximum<double>());
    //maxerr = fmax(maxerr,1.0E-2*tol);
    //cudaErr( cudaPeekAtLastError() );
    //cudaErr( cudaDeviceSynchronize() );
    //printf("Maxerr: %g\n",maxerr);
  
    if (maxerr > tol)
    {
      failcount++;
      // Reject timestep and go back to psi_old
      scale = safety * pow(tol/maxerr,0.25);
      delta_t = scale * delta_t;
      tconst.x = -gamma*delta_t;
      tconst.y = -delta_t;
      //fprintf(stdout,"FAILURE: Time: %g, error: %g, new delta t: %g\n", time, maxerr, delta_t);
    }
    else
    {
      goodcount++;
      t = t + delta_t;
      //delta_t = 0.92 * pow(tol/maxerr,0.20) * delta_t;
      scale = safety * pow(tol/maxerr,0.20);
      if (scale < 2.0)
      {
        delta_t = scale * delta_t;
      }
      else
      {
        delta_t = 2.0 * delta_t;
      }
      tconst.x = -gamma*delta_t;
      tconst.y = -delta_t;
      //fprintf(stdout,"SUCCESS: Time: %g, error: %g, new delta t: %g\n", time, maxerr, delta_t);
      //printf("New dt: %g\n",delta_t);

      // Copy psi to psi_old using copy_kernel
      #ifdef TIMING_MUCH
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);
      #endif
      copyKernel<<< num_blocks, threads_per_block >>>( n, arrs_host[8], arrs_host[0] );
      #ifdef TIMING_MUCH
      cudaEventRecord(stop,0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&kerneltime, start, stop);
      //printf("Elapsed time for copying psi to psi_old: %g ms\n",kerneltime);
      #endif
      cudaErr( cudaPeekAtLastError() );
      cudaErr( cudaDeviceSynchronize() );

      // Update maximum mod of wavefunction
      oldmax = thrust::reduce(rho_thrust, rho_thrust + n, 0.0,thrust::maximum<double>());
      //printf("Oldmax: %g\n",oldmax);
  
      // Check for exit
      if (can_has_exit)
      {
         delta_t = saved_delta_t;
         break;
      }
  
    }
  
  }
  //printf("Good: %d, Bad: %d\n",goodcount,failcount);

  trans_err = cufftExecZ2Z(plan, arrs_host[0], arrs_host[0], CUFFT_INVERSE);
  if (trans_err != CUFFT_SUCCESS)
  {
    printf("Error on final transform to position %d: %d",i,trans_err);
    return;
  }
  #ifdef TIMING_MUCH
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  #endif
  mulrcKernel<<< num_blocks, threads_per_block >>>( n, renorm_pspace, arrs_host[0] );
  #ifdef TIMING_MUCH
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kerneltime, start, stop);
  printf("Elapsed time for renormalizing psi[0]: %g ms\n",kerneltime);
  #endif
  cudaErr( cudaPeekAtLastError() );
  cudaErr( cudaDeviceSynchronize() );

  #ifdef TIMING_MUCH
  for (ii=0; ii<5; ii++)
  {
    for (jj=0; jj<7; jj++)
    {
      timearray[jj][ii] = timearray[jj][ii] / 100.0;
    }
  }
  for (jj=0; jj<6; jj++)
  {
    printf("Average time for stage %d transform to position: %g ms\n", jj+1, timearray[jj][0]);
    printf("Average time for stage %d position space kernel: %g ms\n", jj+1, timearray[jj][1]);
    printf("Average time for stage %d transform to momentum: %g ms\n", jj+1, timearray[jj][2]);
    printf("Average time for stage %d momentum space kernel: %g ms\n", jj+1, timearray[jj][3]);
    printf("Average time for stage %d memory shuffle kernel: %g ms\n", jj+1, timearray[jj][4]);
  }
  printf("Average time for abs kernel: %g ms\n", timearray[6][0]);
  printf("Average time for error kernel: %g ms\n", timearray[6][1]);
  #endif
  
}


