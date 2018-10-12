include ../Makeinit

LIBPCGCOARRAYROOT=../lib/
LIBPCGCOARRAY=libcoarraysolver.a

COARRAYENABLE=#shared

F90FLAGSUB=
coarraytype=
toutput=0
ifeq ($(COARRAYENABLE),dist)
 coarraytype=.dist
 F90FLAGS+=-I$(LIBPCGCOARRAYROOT)
 F90FLAGSUB=-coarray=distributed
 toutput=1
endif

ifeq ($(COARRAYENABLE),shared)
 coarraytype=.shared
 F90FLAGS+=-I$(LIBPCGCOARRAYROOT)
 F90FLAGSUB=-coarray=shared -coarray-num-images=4
 toutput=2
endif

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
 debugext=_debug
else
 DEBUG=
 debugext=
endif


EXEC=prepafile

ifneq ($(COARRAYENABLE),)
 EXEC=pcgrowcoarray
endif

all: $(EXEC)

prepafile: prepafile.f90 
	$(F90) $(F90FLAGS) $^ -o $@$(debugext) $(LIBSPARSEROOT)$(LIBSPARSE)

pcgrowcoarray: pcgrowcoarray.f90
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) -fpp -DTOUTPUT=$(toutput) $^ -o $@$(coarraytype)$(debugext) $(LIBSPARSEROOT)$(LIBSPARSE) $(LIBPCGCOARRAYROOT)$(LIBPCGCOARRAY)

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod *.dist* *.shared* prepafile

