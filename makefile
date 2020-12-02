OBJFILES = getCard.o absArrayKernel.o mulrcKernel.o copyKernel.o errorKernel.o initialProjectorKernel.o pSpaceKernel.o kSpaceKernel.o shuffleKernel.o HDF5io.o PGPEvolver.o main.o

program : describeVersion.h $(OBJFILES)
	nvcc -O3 -I. -I/usr/include -arch=sm_35 --ptxas-options=-v -maxrregcount=32 $(OBJFILES) -L/usr/lib64 -lcuda -lcudart -lcufft -lhdf5 -lz
	#rm describeVersion.h

%.o : %.cu
	nvcc -c -O3 -I.  -arch=sm_35 --ptxas-options=-v -maxrregcount=32 $< -o $@

%.o : %.cpp
	nvcc -c -O3  -I.  -arch=sm_35 --ptxas-options=-v -maxrregcount=32 $< -o $@

describeVersion.h :
	./describeVersion.sh

clean :
	rm *.o *.out

