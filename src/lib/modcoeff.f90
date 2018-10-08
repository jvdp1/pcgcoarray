module modcoeff
 !$ use omp_lib
 use modkind
 !use modsparse
 implicit none
 private
 public::gen_coeff

 !ABSTRACTS
 type,abstract::gen_coeff
  contains
  private
  procedure(multbyv_gen),public,deferred::multbyv
 end type

 abstract interface
  subroutine multbyv_gen(this,x,y,starteq,endeq)
   import::gen_coeff,int4,real8
   class(gen_coeff),intent(inout)::this
   integer(kind=int4),intent(in)::starteq,endeq
   real(kind=real8),intent(in)::x(:)
   real(kind=real8),intent(inout)::y(:)
  end subroutine
 end interface

contains

!PUBLIC

!PRIVATE

end module
