module modsparse
 use modkind
 implicit none
 private
 public::coo,csr

 type::sparse
  private
  integer(kind=int4)::unlog=6
  integer(kind=int4)::n,m      !size of the matrix
  integer(kind=intnel),allocatable,public::ia(:)   !should be private
  integer(kind=int4),allocatable,public::ja(:)     !should be private
  real(kind=real8),allocatable,public::a(:)        !should be private
  contains
   procedure::setoutputunit
 end type

 type,extends(sparse)::coo
  private
  integer(kind=int4)::nel
 contains
  procedure::alloc=>alloc_coo
 end type

 type,extends(sparse)::csr
  private
  integer(kind=int4)::sizeia
 contains
  procedure::get=>getfromfile_csr  
  procedure::alloc=>alloc_csr  
  procedure::get_dimension_1
  procedure::get_dimension_2
  procedure::diagcol=>diagcol_csr  
  procedure::diagrow=>diagrow_csr  
  procedure::sub=>subsparse_csr  
  procedure::multv=>multgenv_csr
  procedure::print=>print_csr
  procedure::printfile=>printtofile_csr
  procedure::printbin=>printtobin_csr
  procedure::reset=>reset_csr
 end type


contains

!SPARSE
subroutine setoutputunit(sparsee,un)
 class(sparse),intent(inout)::sparsee
 integer(kind=int4)::un

 sparsee%unlog=un

end subroutine

!COO
subroutine alloc_coo(sparse,nel,n,m,unlog)
 class(coo),intent(inout)::sparse
 integer(kind=intnel),intent(in)::nel
 integer(kind=int4),intent(in)::n
 integer(kind=int4),intent(in),optional::m
 integer(kind=int4),intent(in),optional::unlog

 sparse%nel=nel
 sparse%n=n
 sparse%m=n
 if(present(m))sparse%m=m
 allocate(sparse%ia(nel),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ja=0
 sparse%a=0.d0
 if(present(unlog))call sparse%setoutputunit(unlog)
end subroutine
 
!CSR
!csr: get matrix
subroutine getfromfile_csr(sparse,namefile,unlog)
 class(csr),intent(inout)::sparse
 integer(kind=int4),intent(in),optional::unlog
 character(len=*),intent(in)::namefile

 integer(kind=int4)::i,j,nel,un
 
 if(present(unlog))call sparse%setoutputunit(unlog)
 write(sparse%unlog,'(/2a)')' Read the matrix from file ',trim(namefile)

 open(newunit=un,file=trim(namefile),status='old',action='read',access='stream')
 read(un)i,j,nel
 call sparse%alloc(nel,i,j)
 
 read(un)sparse%ia
 read(un)sparse%ja
 read(un)sparse%a
 close(un)

 write(sparse%unlog,'(2a/)')' End of reading the matrix from ',trim(namefile)

end subroutine

subroutine alloc_csr(sparse,nel,n,m,unlog)
 class(csr),intent(inout)::sparse
 integer(kind=intnel),intent(in)::nel
 integer(kind=int4),intent(in)::n
 integer(kind=int4),intent(in),optional::m,unlog

 sparse%n=n
 sparse%m=n
 if(present(m))sparse%m=m
 sparse%sizeia=n+1
 allocate(sparse%ia(sparse%sizeia),sparse%ja(nel),sparse%a(nel))
 sparse%ia=0
 sparse%ia(sparse%sizeia)=nel+1
 sparse%ja=0
 sparse%a=0.d0
 if(present(unlog))call sparse%setoutputunit(unlog)
 write(sparse%unlog,'(/a,i0,a,i0)')' Size of the sparse matrix          : ',sparse%n, ' x ',sparse%m
 write(sparse%unlog,'(a,i0)')' Number of elements in sparse matrix: ',sparse%ia(sparse%sizeia)-1
end subroutine
 
function subsparse_csr(sparse,startrow,endrow,startcol,endcol)
 type(csr)::subsparse_csr
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startrow,endrow,startcol,endcol
 
 integer(kind=int4)::newn,newm
 integer(kind=int4)::i,j,k,un
 integer(kind=intnel)::nel
 

 write(sparse%unlog,'(/a,i0,a,i0)')' Extraction of the rows from ',startrow,' to ',endrow
 write(sparse%unlog,'(a,i0,a,i0/)')' Extraction of the columns from ',startcol,' to ',endcol

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
  do j=sparse%ia(i),sparse%ia(i+1)-1
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

!csr: get details
function get_dimension_1(sparse)
 class(csr),intent(in)::sparse
 integer(kind=int4)::get_dimension_1

 get_dimension_1=sparse%n
end function

function get_dimension_2(sparse)
 class(csr),intent(in)::sparse
 integer(kind=int4)::get_dimension_2

 get_dimension_2=sparse%m
end function

!csr: multiplications
subroutine multgenv_csr(sparse,x,y)
 !Computes y=0.d0*y+sparse*x
 class(csr),intent(in)::sparse
 real(kind=real8),intent(in)::x(:)
 real(kind=real8),intent(out)::y(:)

 character(len=1)::matdescra(6)

 matdescra(1)='G'
 matdescra(4)='F'
 
 call mkl_dcsrmv('N',sparse%n,sparse%m,1.d0,matdescra,sparse%a,sparse%ja,sparse%ia(1:sparse%n),sparse%ia(2:sparse%n+1),x,0.d0,y)
 
end subroutine

!csr: print
subroutine print_csr(sparse)
 class(csr),intent(in)::sparse

 integer(kind=int4)::i
 integer(kind=intnel)::j

 do i=1,sparse%n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   write(sparse%unlog,'(i10,i10,f16.8)'),i,sparse%ja(j),sparse%a(j)
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
   write(un,'(i13,i13,f16.8)')i,sparse%ja(j),sparse%a(j)
  enddo
 enddo

end subroutine

subroutine printtobin_csr(sparse,un,typesp)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::un
 character(len=*),intent(in)::typesp

 integer(kind=int4)::i,unin
 integer(kind=intnel)::j
 character(len=10)::cdummy
 
 write(cdummy,'(i0)')un
 open(newunit=unin,file='subpcg.'//adjustl(typesp(:len_trim(typesp)))//adjustl(cdummy(:len_trim(cdummy))),action='write',access='stream',buffered='yes')
 write(unin)sparse%n,sparse%m,sparse%ia(sparse%n+1)-1
 write(unin)sparse%ia
 write(unin)sparse%ja
 write(unin)sparse%a
 close(unin)

end subroutine

!csr: reset
subroutine reset_csr(sparse)
 class(csr),intent(inout)::sparse
 
 sparse%n=0
 sparse%m=0
 sparse%sizeia=0
 deallocate(sparse%ia,sparse%ja,sparse%a)

end subroutine

function diagcol_csr(sparse,startcol)
 real(kind=real8),allocatable::diagcol_csr(:)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startcol
 
 integer(kind=int4)::newn,newm
 integer(kind=int4)::i,j,k,un
 
 write(sparse%unlog,'(/a,i0)')' Extraction of the diagonal elements from ',startcol

 allocate(diagcol_csr(sparse%m))
 diagcol_csr=0.d0
 newn=0
 do i=startcol,startcol+sparse%m-1
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(k+startcol-1.eq.i)then
    newn=newn+1
    diagcol_csr(newn)=sparse%a(j)
   endif
  enddo
 enddo

end function

function diagrow_csr(sparse,startrow)
 real(kind=real8),allocatable::diagrow_csr(:)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startrow
 
 integer(kind=int4)::newn,newm
 integer(kind=int4)::i,j,k,un
 
 write(sparse%unlog,'(/a,i0)')' Extraction of the diagonal elements from ',startrow

 allocate(diagrow_csr(sparse%n))
 diagrow_csr=0.d0
 newn=0
 do i=1,sparse%n
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(k.eq.i+startrow-1)then
    newn=newn+1
    diagrow_csr(newn)=sparse%a(j)
   endif
  enddo
 enddo

end function
end module
