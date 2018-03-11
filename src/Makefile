F90=ifort

F90FLAGS=-O3 -heap-arrays -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -parallel -lpthread

F90FLAGSUB=
coarray=

ifeq ($(COARRAYENABLE),1)
 coarray=coarray
 F90FLAGSUB=-coarray=distributed -coarray-num-images=3
endif

DEBUG=

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
endif

LDFLAGS=

EXEC=prepafile pcg

ifeq ($(COARRAYENABLE),1)
 EXEC=pcg$(coarray) pcgrow$(coarray)
endif

all: $(EXEC)

prepafile: modkind.o modsparse.o prepafile.f90 
	$(F90) $(F90FLAGS) $(LDFLAGS) $^ -o $@

pcg$(coarray): modkind.o modsparse.o pcg$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@

pcgrow$(coarray): modkind.o modsparse.o pcgrow$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) $< -o $@

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod pcg pcgcoarray pcgrowcoarray prepafile

