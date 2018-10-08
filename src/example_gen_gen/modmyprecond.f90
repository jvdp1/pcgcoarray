module modmyprecond
 !$ use omp_lib
 use modkind
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
  integer(kind=int4)::dim1
  real(kind=real8),allocatable::array(:)
 contains
  procedure,public::solve
 end type

contains

!SOLVE
subroutine solve(this,x,y)
 !solve this*x=y, by doing x=inv(this)*y
 class(arrayprecond),intent(inout)::this
 real(kind=real8),intent(out)::x(:)
 real(kind=real8),intent(inout)::y(:)
 
 integer(kind=int4)::i
 logical,save::lnotinverse=.true.

 if(lnotinverse)then
  do i=1,this%dim1
   if(this%array(i).ne.0_real8)this%array(i)=1._real8/this%array(i)
  enddo
  lnotinverse=.false.
 endif
 
 do i=1,this%dim1
  x(i)=this%array(i)*y(i)
 enddo

end subroutine

end module
