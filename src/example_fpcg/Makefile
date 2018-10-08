include ../Makeinit

LIBPCGCOARRAYROOT=../lib/
LIBPCGCOARRAY=libpcgcoarray.a

COARRAYENABLE=shared

F90FLAGSUB=

coarraytype=
toutput=0
coarray=0
ifeq ($(COARRAYENABLE),dist)
 coarraytype=.dist
 F90FLAGS+=-I$(LIBPCGCOARRAYROOT)
 F90FLAGSUB=-coarray=distributed
 toutput=1
 coarray=1
endif

ifeq ($(COARRAYENABLE),shared)
 coarraytype=.shared
 F90FLAGS+=-I$(LIBPCGCOARRAYROOT)
 F90FLAGSUB=-coarray=shared -coarray-num-images=4
 toutput=2
 coarray=1
endif

ifeq ($(DEBUGENABLE),1)
 DEBUG+=-g -check all -traceback
 debugext=_debug
else
 DEBUG=
 debugext=
endif

OBJ=modmvlr.o
ifneq ($(COARRAYENABLE),)
 OBJ+=modmyprecond.o
endif

ifneq ($(COARRAYENABLE),)
 EXEC=fpcgrowcoarray_mvlr
endif

all: $(EXEC)

fpcgrowcoarray_mvlr: $(OBJ) fpcgrowcoarray_mvlr.f90
	$(F90) $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) -fpp -DTOUTPUT=$(toutput) -DCOARRAY=$(coarray) $^ -o $@$(coarraytype)$(debugext) $(LIBSPARSEROOT)$(LIBSPARSE) $(LIBPCGCOARRAYROOT)$(LIBPCGCOARRAY)

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(F90FLAGSUB) $(DEBUG) -fpp -DTOUTPUT=$(toutput) -DCOARRAY=$(coarray) $< -o $@

cleanlib:
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod *.dist* *.shared* 

