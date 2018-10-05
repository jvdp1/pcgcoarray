# Preconditioned Conjugate Gradient method implemented in Coarray Fortran  


Coarray Fortran library and programs for running the PCG method in parallel.  

The library and programs were compiled and tested with Intel Fortran 2017 and 2018. The Intel MKL library (i.e., BLAS, Sparse BLAS, LAPACK, and Pardiso) was used. Sparse matrices are handled with the library libsparse (https://github.com/jvdp1/libsparse).  

Earlier versions of the programs were compiled with the GNU compiler Collection (GCC) Fortran and OpenCoarrays.  

The subroutine *pcgrowcoarray* performs the PCG method. The coefficient matrix and the preconditioner can be provided to the subroutine *pcgrowcoarray* under different derived types. Different examples are provided. For example, one example uses a coefficient matrix stored in a CRS format, while another example uses a matrix-free coefficient matrix.  


## Compilation  


For compilation of Fortran programs: `make`  
For compilation of shared-memory Coarray Fortran programs: `make COARRAYENABLE=shared`  
For compilation of distributed-memory Coarray Fortran programs: `make COARRAYENABLE=dist`  

Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  


