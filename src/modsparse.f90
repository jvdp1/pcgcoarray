module modsparse
 use modkind
 implicit none
 private
 public::coo,csr

 type::sparse
  integer(kind=int4)::n,m      !size of the matrix
  integer(kind=intnel),allocatable::ia(:)
  integer(kind=int4),allocatable::ja(:)
  real(kind=real8),allocatable::a(:)
 end type

 type,extends(sparse)::coo
  integer(kind=int4)::nel
 contains
  procedure::alloc=>alloc_coo
 end type

 type,extends(sparse)::csr
  integer(kind=int4)::sizeia
 contains
  procedure::alloc=>alloc_csr  
  procedure::diag=>diag_csr  
  procedure::sub=>subsparse_csr  
  procedure::print=>print_csr
  procedure::printfile=>printtofile_csr
  procedure::printbin=>printtobin_csr
  procedure::reset=>reset_csr
 end type


contains

subroutine alloc_coo(sparse,nel,n,m)
 class(coo),intent(inout)::sparse
 integer(kind=intnel),intent(in)::nel
 integer(kind=int4),intent(in)::n
 integer(kind=int4),intent(in),optional::m

 sparse%nel=nel
 sparse%n=n
 sparse%m=n
 if(present(m))sparse%m=m
 allocate(sparse%ia(nel),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ja=0
 sparse%a=0.d0
end subroutine
 
subroutine alloc_csr(sparse,nel,n,m)
 class(csr),intent(inout)::sparse
 integer(kind=intnel),intent(in)::nel
 integer(kind=int4),intent(in)::n
 integer(kind=int4),intent(in),optional::m

 sparse%n=n
 sparse%m=n
 if(present(m))sparse%m=m
 sparse%sizeia=n+1
 allocate(sparse%ia(sparse%sizeia),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%sizeia)=nel+1
 sparse%ja=0
 sparse%a=0.d0
 write(*,'(/a,i0,a,i0)')' Size of the sparse matrix          : ',sparse%n, ' x ',sparse%m
 write(*,'(a,i0)')' Number of elements in sparse matrix: ',sparse%ia(sparse%sizeia)-1
end subroutine
 
function subsparse_csr(sparse,startrow,endrow,startcol,endcol)
 type(csr)::subsparse_csr
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startrow,endrow,startcol,endcol
 
 integer(kind=int4)::newn,newm
 integer(kind=int4)::i,j,k,un
 integer(kind=intnel)::nel
 

 write(*,'(/a,i0,a,i0)')' Extraction of the rows from ',startrow,' to ',endrow
 write(*,'(a,i0,a,i0/)')' Extraction of the columns from ',startcol,' to ',endcol

 newn=endrow-startrow+1
 newm=endcol-startcol+1
 
 nel=0
 do i=startrow,endrow
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(k.ge.startcol.and.k.le.endcol)then
    nel=nel+1
   endif
  enddo
 enddo

 call subsparse_csr%alloc(nel,newn,newm)

 nel=0
 subsparse_csr%ia=0
 
 un=0
 nel=1
 do i=startrow,endrow
  un=un+subsparse_csr%ia(nel)
  nel=nel+1
  do j=sparse%ia(i),sparse%ia(nel)-1
    k=sparse%ja(j)
    if(k.ge.startcol.and.k.le.endcol)then
     subsparse_csr%ia(nel)=subsparse_csr%ia(nel)+1
     subsparse_csr%ja(un+subsparse_csr%ia(nel))=k-startcol+1
     subsparse_csr%a(un+subsparse_csr%ia(nel))=sparse%a(j)
    endif
  enddo
 enddo
 subsparse_csr%ia(1)=1
 do i=2,subsparse_csr%n+1
  subsparse_csr%ia(i)=subsparse_csr%ia(i)+subsparse_csr%ia(i-1)
 enddo

end function

subroutine print_csr(sparse)
 class(csr),intent(in)::sparse

 integer(kind=int4)::i
 integer(kind=intnel)::j

 do i=1,sparse%n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   write(*,'(i10,i10,f16.8)'),i,sparse%ja(j),sparse%a(j)
  enddo
 enddo

end subroutine

subroutine printtofile_csr(sparse,un)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::un

 integer(kind=int4)::i
 integer(kind=intnel)::j

 do i=1,sparse%n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   write(un,'(i10,i10,f16.8)'),i,sparse%ja(j),sparse%a(j)
  enddo
 enddo

end subroutine

subroutine printtobin_csr(sparse,un)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::un

 integer(kind=int4)::i,unin
 integer(kind=intnel)::j
 character(len=10)::cdummy
 
 write(cdummy,'(i0)')un
 open(newunit=unin,file='subpcg.'//adjustl(cdummy(:len_trim(cdummy))),action='write',access='stream',buffered='yes')
 write(unin)sparse%n,sparse%m,sparse%ia(sparse%n+1)-1
 write(unin)sparse%ia
 write(unin)sparse%ja
 write(unin)sparse%a
 close(unin)

end subroutine

subroutine reset_csr(sparse)
 class(csr),intent(inout)::sparse
 
 sparse%n=0
 sparse%m=0
 sparse%sizeia=0
 deallocate(sparse%ia,sparse%ja,sparse%a)

end subroutine

function diag_csr(sparse,startcol)
 real(kind=real8),allocatable::diag_csr(:)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startcol
 
 integer(kind=int4)::newn,newm
 integer(kind=int4)::i,j,k,un
 integer(kind=intnel)::nel
 
 write(*,'(/a,i0)')' Extraction of the diagonal elements from ',startcol

 allocate(diag_csr(sparse%m))
 diag_csr=0.d0
 newn=0
 do i=startcol,startcol+sparse%m-1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(k+startcol-1.eq.i)then
    newn=newn+1
    diag_csr(newn)=sparse%a(j)
   endif
  enddo
 enddo

end function

end module
