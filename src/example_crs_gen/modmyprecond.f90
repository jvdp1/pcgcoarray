module modmyprecond
 use iso_fortran_env
 !$ use omp_lib
#if (COARRAY==1)
 use modprecond,only:gen_precond
#endif
 implicit none
 private
 public::arrayprecond

#if (COARRAY==1)
 type,extends(gen_precond)::arrayprecond
#else
 type::arrayprecond
#endif
  integer(kind=int32)::dim1
  real(kind=real64),allocatable::array(:)
 contains
  procedure,public::solve
 end type

contains

!SOLVE
subroutine solve(this,x,y)
 !solve this*x=y, by doing x=inv(this)*y
 class(arrayprecond),intent(inout)::this
 real(kind=real64),intent(out)::x(:)
 real(kind=real64),intent(inout)::y(:)
 
 integer(kind=int32)::i
 logical,save::lnotinverse=.true.

 if(lnotinverse)then
  do i=1,this%dim1
   if(this%array(i).ne.0_real64)this%array(i)=1._real64/this%array(i)
  enddo
  lnotinverse=.false.
 endif
 
 do i=1,this%dim1
  x(i)=this%array(i)*y(i)
 enddo

end subroutine

end module
