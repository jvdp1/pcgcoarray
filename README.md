# Preconditioned Conjugate Gradient method implemented in Coarray Fortran  


Basic Coarray Fortran programs for running the PCG method in parallel.  

The programs were compiled and tested with Intel Fortran 2017 and 2018. The Intel MKL library (i.e., BLAS, Sparse BLASE, LAPACK, and Pardiso) is used.  
Earlier versions of the programs were compiled with the GNU compiler Collection (GCC) Fortran and OpenCoarrays.  

The program `pcgrowcoarray.f90` requires that the preconditioner is provided as a upper triangular matrix in the CSR format. The system of equations involving the preconditioner is solved using MKL Pardiso.  


## Compilation  


For compilation of Fortran programs: `make`  
For compilation of shared-memory Coarray Fortran programs: `make COARRAYENABLE=shared`  
For compilation of distributed-memory Coarray Fortran programs: `make COARRAYENABLE=dist`  

Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  





