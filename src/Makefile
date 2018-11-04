all: libcoarray example_crs_crs example_crs_gen example_fpcg example_gen_gen

lib: libcoarray

libcoarray:	
	cd lib;make

example_crs_crs: lib
	cd example_crs_crs;make;make cleanlib;make COARRAYENABLE=shared

example_crs_gen: lib
	cd  example_crs_gen;make;make cleanlib;make COARRAYENABLE=shared

example_fpcg: lib
	cd example_fpcg;make;make cleanlib;make COARRAYENABLE=shared

example_gen_gen: lib
	cd example_gen_gen;make;make cleanlib;make COARRAYENABLE=shared

clean:
	cd lib;make clean
	cd example_crs_crs;make clean
	cd example_crs_gen;make clean
	cd example_fpcg;make clean
	cd example_gen_gen;make clean

cleanlib:
	rm -rf */*.o */*.mod
