F90=ifort

F90FLAGS=-O3 -heap-arrays -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -parallel -lpthread

F90FLAGSUB=
coarray=

ifeq ($(COARRAYENABLE),1)
 coarray=coarray
 F90FLAGSUB=-coarray=distributed -coarray-num-images=3
endif

LDFLAGS=


EXEC=prepafile pcg$(coarray)

all: $(EXEC)

prepafile: modkind.o modsparse.o prepafile.f90 
	$(F90) $(F90FLAGS) $(LDFLAGS) $^ -o $@

pcg$(coarray): modkind.o modsparse.o pcg$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $^ -o $@

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $< -o $@

.PHONY: proper

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod pcg pcgcoarray prepafile

