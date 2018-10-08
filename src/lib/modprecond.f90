module modprecond
 !$ use omp_lib
 use modkind
 !use modsparse
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
   import::gen_precond,real8
   class(gen_precond),intent(inout)::this
   real(kind=real8),intent(out)::x(:)
   real(kind=real8),intent(inout)::y(:)
  end subroutine
 end interface

contains

!PUBLIC

!PRIVATE

end module
