module modcoeff
 use iso_fortran_env
 !$ use omp_lib
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
   import::gen_coeff,int32,real64
   class(gen_coeff),intent(inout)::this
   integer(kind=int32),intent(in)::starteq,endeq
   real(kind=real64),intent(in)::x(:)
   real(kind=real64),intent(inout)::y(:)
  end subroutine
 end interface

contains

!PUBLIC

!PRIVATE

end module
