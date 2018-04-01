F90=ifort

F90FLAGS=-O3 -heap-arrays -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -parallel -lpthread

F90FLAGSUB=
coarraytype=
coarray=

ifeq ($(COARRAYENABLE),dist)
 coarray=coarray
 coarraytype=.dist
 F90FLAGSUB=-coarray=distributed
endif

ifeq ($(COARRAYENABLE),shared)
 coarray=coarray
 coarraytype=.shared
 F90FLAGSUB=-coarray=shared -coarray-num-images=4
endif

DEBUG=
debugext=

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
 debugext=_debug
endif

LDFLAGS=

EXEC=prepafile pcg

ifneq ($(COARRAYENABLE),)
 EXEC=pcg$(coarray) pcgrow$(coarray)
endif

all: $(EXEC)

prepafile: modkind.o modsparse.o prepafile.f90 
	$(F90) $(F90FLAGS) $(LDFLAGS) $^ -o $@$(debugext)

pcg$(coarray): modkind.o modsparse.o pcg$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@$(coarraytype)$(debugext)

pcgrow$(coarray): modkind.o modsparse.o pcgrow$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@$(coarraytype)$(debugext)

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) $< -o $@

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod pcg *.dist* *.shared* prepafile

