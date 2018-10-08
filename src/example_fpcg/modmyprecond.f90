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
  real(kind=real8),allocatable::array1(:)
  real(kind=real8),allocatable::array2(:)
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
 integer(kind=int4)::iter=1
 logical,save::lnotinverse1=.true.
 logical,save::lnotinverse2=.true.


 if(iter.eq.1)then
  if(lnotinverse1)then
   call inversediag(this%dim1,this%array1)
   lnotinverse1=.false.
  endif
  call solvediag(this%dim1,this%array1,x,y)
  iter=2
 elseif(iter.eq.2)then
  if(lnotinverse2)then
   call inversediag(this%dim1,this%array2)
   lnotinverse2=.false.
  endif
  call solvediag(this%dim1,this%array2,x,y)
  iter=1
 endif

end subroutine

subroutine inversediag(dim1,array)
 integer(kind=int4),intent(in)::dim1
 real(kind=real8),intent(inout)::array(:)

 integer(kind=int4)::i

 do i=1,dim1
  if(array(i).ne.0_real8)array(i)=1._real8/array(i)
 enddo

end subroutine

subroutine solvediag(dim1,array,x,y)
 integer(kind=int4),intent(in)::dim1
 real(kind=real8),intent(inout)::array(:)
 real(kind=real8),intent(inout)::y(:)
 real(kind=real8),intent(out)::x(:)

 integer(kind=int4)::i

 do i=1,dim1
  x(i)=array(i)*y(i)
 enddo

end subroutine

end module
