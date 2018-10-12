module modcoarrayop
 use iso_fortran_env
 implicit none
 private
 public::co_sum

 interface co_sum
  module procedure co_sum_scal_real64
 end interface

contains

subroutine co_sum_scal_real64(this)
 real(kind=real64),intent(inout)::this[*]

 integer(kind=int32)::i

 sync all

 if(this_image().eq.1)then
  !sum to image 1
  do i=2,num_images()
   this=this+this[i]
  enddo
 endif

 sync all
 
 !update all other images
 if(this_image().ne.1)this=this[1]

 sync all

end subroutine

end module

