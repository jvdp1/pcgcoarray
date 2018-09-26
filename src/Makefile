F90=ifort

LIBSPARSEROOT=../../
LIBSPARSE=libsparse/src/lib/

F90FLAGS=-O3 -heap-arrays -qopenmp -parallel -stand f08
F90FLAGS+=-I$(LIBSPARSEROOT)$(LIBSPARSE)
F90FLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64
F90FLAGS+=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
F90FLAGS+=-lpthread -ldl

F90FLAGSUB=
coarraytype=

toutput=0

PARDISOENABLE=1

ifeq ($(PARDISOENABLE),1)
 PARDISO=1
 F90FLAGS+=
endif

ifeq ($(COARRAYENABLE),dist)
 coarraytype=.dist
 F90FLAGSUB=-coarray=distributed
 toutput=1
endif

ifeq ($(COARRAYENABLE),shared)
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
 EXEC=pcgrowcoarray
endif

all: $(EXEC)

prepafile: modkind.o $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a prepafile.f90 
	$(F90) $(F90FLAGS) $(LDFLAGS) $^ -o $@$(debugext) $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a

pcgrowcoarray: modkind.o pcgrowcoarray.o
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(LDFLAGS) $(DEBUG) $^ -o $@$(coarraytype)$(debugext) $(LIBSPARSEROOT)$(LIBSPARSE)/libsparse.a

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) -fpp -DTOUTPUT=$(toutput) -D_PARDISO=$(PARDISO) $< -o $@

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod pcg *.dist* *.shared* prepafile

