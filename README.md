# Preconditioned Conjugate Gradient method implemented in Coarray Fortran  


Basic Coarray Fortran programs for running the PCG method in parallel.  

For compilation of Fortran programs: `make`  
For compilation of shared memory coarray Fortran programs: `make COARRAYENABLE=shared`  
For compilation of distributed memory coarray Fortran programs: `make COARRAYENABLE=dist`  

Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  





