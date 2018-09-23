F90=ifort

F90FLAGS=-O3 -heap-arrays -qopenmp -parallel -stand f08
F90FLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64
F90FLAGS+=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
F90FLAGS+=-lpthread -ldl

F90FLAGSUB=
coarraytype=
coarray=

toutput=0

PARDISOENABLE=1

ifeq ($(PARDISOENABLE),1)
 PARDISO=1
 F90FLAGS+=
endif

ifeq ($(COARRAYENABLE),dist)
 coarray=coarray
 coarraytype=.dist
 F90FLAGSUB=-coarray=distributed
 toutput=1
endif

ifeq ($(COARRAYENABLE),shared)
 coarray=coarray
 coarraytype=.shared
 F90FLAGSUB=-coarray=shared -coarray-num-images=4
 toutput=2
endif

DEBUG=
debugext=

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
 debugext=_debug
endif

LDFLAGS=

EXEC=prepafile

ifneq ($(COARRAYENABLE),)
 EXEC=pcgrow$(coarray)
endif

all: $(EXEC)

prepafile: modkind.o modsparse_old.o prepafile.f90 
	$(F90) $(F90FLAGS) $(LDFLAGS) $^ -o $@$(debugext)

pcgrow$(coarray): modkind.o modsparse_old.o pcgrow$(coarray).o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@$(coarraytype)$(debugext)

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) -fpp -DTOUTPUT=$(toutput) -D_PARDISO=$(PARDISO) $< -o $@

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod pcg *.dist* *.shared* prepafile

