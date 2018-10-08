include ../Makeinit

COARRAYENABLE=shared

coarray=0
coarraytype=
toutput=0
ifeq ($(COARRAYENABLE),dist)
 coarraytype=.dist
 F90FLAGSUB=-coarray=distributed
 coarray=1
 toutput=1
endif

ifeq ($(COARRAYENABLE),shared)
 coarraytype=.shared
 F90FLAGSUB=-coarray=shared -coarray-num-images=4
 coarray=1
 toutput=2
endif

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
 debugext=_debug
else
 DEBUG=
 debugext=
endif


all:	libpcgcoarray.a

OBJ=modkind.o modprecond.o modcoeff.o modcoarraysolver.o 

libpcgcoarray.a: $(OBJ)
	ar cr  libpcgcoarray.a $(OBJ)

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) $< -o $@

cleanlib:
	rm -rf *.o *.mod *.a

clean:
	rm -rf *.o *.mod *.a *.dist* *.shared* 

