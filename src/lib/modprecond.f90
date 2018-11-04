module modprecond
 use iso_fortran_env
 !$ use omp_lib
 implicit none
 private
 public::gen_precond

 !ABSTRACTS
 type,abstract::gen_precond
  contains
  private
  procedure(solve_gen),public,deferred::solve
 end type

 abstract interface
  subroutine solve_gen(this,x,y)
   import::gen_precond,real64
   class(gen_precond),intent(inout)::this
   real(kind=real64),intent(out)::x(:)
   real(kind=real64),intent(inout)::y(:)
  end subroutine
 end interface

contains

!PUBLIC

!PRIVATE

end module
