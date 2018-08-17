#if (_PARDISO==1)
include 'mkl_pardiso.f90'
#endif

module modsparse
 use modkind
#if (_PARDISO==1)
 use mkl_pardiso
#endif
 !$ use omp_lib
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
  logical::lsquare
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
  procedure::sort=>sort_csr  
  procedure::diagcol=>diagcol_csr  
  procedure::diagrow=>diagrow_csr  
  procedure::sub=>subsparse_csr  
  procedure::subup=>subsparse_up_csr  
  procedure::subdiag=>subsparse_diag_csr  
  procedure::subtriu=>subsparse_triu_csr  
  procedure::multv=>multgenv_csr
  procedure::solve=>solve_csr
  procedure::print=>print_csr
  procedure::printfile=>printtofile_csr
  procedure::printbin=>printtobin_csr
  procedure::reset=>reset_csr     !possibility to put it as a final procedure?
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
 sparse%lsquare=.true.
 if(sparse%n.ne.sparse%m)sparse%lsquare=.false.
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
 sparse%lsquare=.true.
 if(sparse%n.ne.sparse%m)sparse%lsquare=.false.
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

function subsparse_triu_csr(sparse,startrow,endrow,startcol,endcol,oaway)
 type(csr)::subsparse_triu_csr
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::startrow,endrow,startcol,endcol
 integer(kind=int4),intent(in),optional::oaway
 
 integer(kind=int4)::newn,newm,away
 integer(kind=int4)::i,j,k,un
 integer(kind=intnel)::nel
 logical::lidentity
 

 write(sparse%unlog,'(/a,i0,a,i0)')' Extraction of the rows from ',startrow,' to ',endrow
 write(sparse%unlog,'(a,i0,a,i0)')' Extraction of the columns from ',startcol,' to ',endcol

 newn=endrow-startrow+1
 newm=endcol-startcol+1

 lidentity=.false.
 away=max(endrow,endcol)
 if(present(oaway))away=oaway 
 if(away<0)then
  away=0
  lidentity=.true.
 endif

 write(sparse%unlog,'(a,i0)')' Number of diagonals above the main diagonal: ',away
 
 nel=0
 do i=startrow,endrow
  do j=sparse%ia(i),sparse%ia(i+1)-1
   k=sparse%ja(j)
   if(k.ge.startcol.and.k.le.endcol.and.k.ge.i.and.k.le.i+away)then
    nel=nel+1
   endif
  enddo
 enddo

 call subsparse_triu_csr%alloc(nel,newn,newm)

 nel=0
 subsparse_triu_csr%ia=0
 
 un=0
 nel=1
 do i=startrow,endrow
  un=un+subsparse_triu_csr%ia(nel)
  nel=nel+1
  do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startcol.and.k.le.endcol.and.k.ge.i.and.k.le.i+away)then
     subsparse_triu_csr%ia(nel)=subsparse_triu_csr%ia(nel)+1
     subsparse_triu_csr%ja(un+subsparse_triu_csr%ia(nel))=k-startcol+1
     subsparse_triu_csr%a(un+subsparse_triu_csr%ia(nel))=sparse%a(j)
    endif
  enddo
 enddo
 if(lidentity)subsparse_triu_csr%a=1.d0
 subsparse_triu_csr%ia(1)=1
 do i=2,subsparse_triu_csr%n+1
  subsparse_triu_csr%ia(i)=subsparse_triu_csr%ia(i)+subsparse_triu_csr%ia(i-1)
 enddo

end function

function subsparse_up_csr(sparse,startrow,endrow,startcol,endcol)
 type(csr)::subsparse_up_csr
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
   if(k.ge.startcol.and.k.le.endcol.and.k.ge.i)then
   !if(k.ge.startcol.and.k.le.endcol.and.k.ge.i.and.k.lt.i+100)then
    nel=nel+1
   endif
  enddo
 enddo

 call subsparse_up_csr%alloc(nel,newn,newm)

 nel=0
 subsparse_up_csr%ia=0
 
 un=0
 nel=1
 do i=startrow,endrow
  un=un+subsparse_up_csr%ia(nel)
  nel=nel+1
  do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startcol.and.k.le.endcol.and.k.ge.i)then
    !if(k.ge.startcol.and.k.le.endcol.and.k.ge.i.and.k.lt.i+100)then
     subsparse_up_csr%ia(nel)=subsparse_up_csr%ia(nel)+1
     subsparse_up_csr%ja(un+subsparse_up_csr%ia(nel))=k-startcol+1
     subsparse_up_csr%a(un+subsparse_up_csr%ia(nel))=sparse%a(j)
    endif
  enddo
 enddo
 subsparse_up_csr%ia(1)=1
 do i=2,subsparse_up_csr%n+1
  subsparse_up_csr%ia(i)=subsparse_up_csr%ia(i)+subsparse_up_csr%ia(i-1)
 enddo

end function

function subsparse_diag_csr(sparse,startrow,endrow,startcol,endcol)
 type(csr)::subsparse_diag_csr
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
   if(k.ge.startcol.and.k.le.endcol.and.k.eq.i)then
    nel=nel+1
   endif
  enddo
 enddo

 call subsparse_diag_csr%alloc(nel,newn,newm)

 nel=0
 subsparse_diag_csr%ia=0
 
 un=0
 nel=1
 do i=startrow,endrow
  un=un+subsparse_diag_csr%ia(nel)
  nel=nel+1
  do j=sparse%ia(i),sparse%ia(i+1)-1
    k=sparse%ja(j)
    if(k.ge.startcol.and.k.le.endcol.and.k.eq.i)then
     subsparse_diag_csr%ia(nel)=subsparse_diag_csr%ia(nel)+1
     subsparse_diag_csr%ja(un+subsparse_diag_csr%ia(nel))=k-startcol+1
     subsparse_diag_csr%a(un+subsparse_diag_csr%ia(nel))=sparse%a(j)
    endif
  enddo
 enddo
 subsparse_diag_csr%ia(1)=1
 do i=2,subsparse_diag_csr%n+1
  subsparse_diag_csr%ia(i)=subsparse_diag_csr%ia(i)+subsparse_diag_csr%ia(i-1)
 enddo

end function

subroutine sort_csr(sparse)
 ! sort vectors ja and a by increasing order
 class(csr),intent(inout)::sparse

 integer(kind=int4)::dir,endd,i,j,k,n,start,stkpnt
 integer(kind=int4)::d1,d2,d3,dmnmx,tmp
 integer(kind=int4)::stack(2,32)
 integer(kind=int4),allocatable::d(:)
 integer(kind=int4),parameter::select=20
 real(kind=real8)::umnmx,tmpu
 real(kind=real8),allocatable::u(:)

 do k=1,sparse%n
  n=sparse%ia(k+1)-sparse%ia(k)
  if( n > 1 ) then
   allocate(d(n),u(n))
   !copy of the vector to be sorted
   d=sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)
   u=sparse%a(sparse%ia(k):sparse%ia(k+1)-1)
   !sort the vectors
   !from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Quick return if possible
   stkpnt = 1
   stack( 1, 1 ) = 1
   stack( 2, 1 ) = n
   10 start = stack( 1, stkpnt )
   endd = stack( 2, stkpnt )
   stkpnt = stkpnt - 1
   IF( endd-start <= select .AND. endd-start > 0 ) THEN
   !Do Insertion sort on D( START:ENDD )
   !Sort into increasing order
     DO i = start + 1, endd
      DO j = i, start + 1, -1
       IF( d( j ) < d( j-1 ) ) THEN
         dmnmx = d( j )
         d( j ) = d( j-1 )
         d( j-1 ) = dmnmx
         umnmx = u( j )
         u( j ) = u( j-1 )
         u( j-1 ) = umnmx
       ELSE
         CYCLE
       END IF
      END DO
     END DO
   ELSE IF( endd-start > select ) THEN
   !Partition D( START:ENDD ) and stack parts, largest one first
   !Choose partition entry as median of 3
    d1 = d( start )
    d2 = d( endd )
    i = ( start + endd ) / 2
    d3 = d( i )
    IF( d1 < d2 ) THEN
      IF( d3 < d1 ) THEN
        dmnmx = d1
      ELSE IF( d3 < d2 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d2
      END IF
    ELSE
      IF( d3 < d2 ) THEN
        dmnmx = d2
      ELSE IF( d3 < d1 ) THEN
        dmnmx = d3
      ELSE
        dmnmx = d1
      END IF
    END IF
    !Sort into increasing order
     i = start - 1
     j = endd + 1
     90 j = j - 1
     IF( d( j ) > dmnmx ) GO TO 90
     110 i = i + 1
     IF( d( i ) < dmnmx ) GO TO 110
     IF( i < j ) THEN
       tmp = d( i )
       d( i ) = d( j )
       d( j ) = tmp
       tmpu = u( i )
       u( i ) = u( j )
       u( j ) = tmpu
       GO TO 90
     END IF
     IF( j-start > endd-j-1 ) THEN
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
     ELSE
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = j + 1
       stack( 2, stkpnt ) = endd
       stkpnt = stkpnt + 1
       stack( 1, stkpnt ) = start
       stack( 2, stkpnt ) = j
     END IF
   END IF
   IF( stkpnt > 0 ) GO TO 10
   !end from dlasrt.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !copy back the sorted vectors
   sparse%ja(sparse%ia(k):sparse%ia(k+1)-1)=d
   sparse%a(sparse%ia(k):sparse%ia(k+1)-1)=u
   deallocate(d,u)
  endif
 enddo
 
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

!csr: multiplications and solve
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

#if (_PARDISO==1)
subroutine solve_csr(sparse,x,y)
 !sparse*x=y
 class(csr),intent(in)::sparse
 real(kind=real8),intent(out)::x(:)
 real(kind=real8),intent(inout)::y(:)

 !Pardiso variables
 integer(kind=int4)::mtype=-2
 !integer(kind=int4)::mtype=11
 integer(kind=int4)::solver=0,error,phase,maxfct=1,mnum=1,nrhs=1
 integer(kind=int4)::idum(1)
 integer(kind=int4),save::iparm(64),msglvl=1
 real(kind=real8)::ddum(1)
 type(MKL_PARDISO_HANDLE),allocatable,save::pt(:)
 logical,save::lpardisofirst=.true.

 integer(kind=int4)::i
 !$ real(kind=real8)::t1

 if(.not.sparse%lsquare)then
  write(sparse%unlog,'(a)')' Warning: the sparse matrix is not squared!'
  return
 endif

 if(lpardisofirst)then
  !$ t1=omp_get_wtime()
  !Preparation of Cholesky of A11 with Pardiso
  !initialize pt
  allocate(pt(64))
  do i = 1, 64
    pt(i)%DUMMY=0 
  enddo
 
  !initialize iparm
  call pardisoinit(pt,mtype,iparm)
!  do i=1,64
!   write(sparse%unlog,*)'iparm',i,iparm(i)
!  enddo
  
  !Ordering and factorization
  phase=12
  iparm(2)=3
  iparm(27)=1
  write(sparse%unlog,'(a)')' Start ordering and factorization'
  call pardiso(pt,maxfct,mnum,mtype,phase,&
               sparse%get_dimension_1(),sparse%a,sparse%ia,sparse%ja,&
               idum,nrhs,iparm,msglvl,ddum,ddum,error)
  call checkparido(phase,error) 
 
  write(sparse%unlog,'(a,i0)')' Number of nonzeros in factors  = ',iparm(18)
  write(sparse%unlog,'(a,i0)')' Number of factorization MFLOPS = ',iparm(19)
  !$ write(sparse%unlog,'(a,i0)')' Elapsed time                   = ',omp_get_wtime()-t1
 endif 

 !Solving
 phase=33
 iparm(27)=0
 call pardiso(pt,maxfct,mnum,mtype,phase,&
              sparse%get_dimension_1(),sparse%a,sparse%ia,sparse%ja,&
              idum,nrhs,iparm,msglvl,y,x,error)
 call checkparido(phase,error) 

 msglvl=0
 lpardisofirst=.false.

contains

 subroutine checkparido(phase,error)
  integer(kind=int4),intent(in)::phase,error
  if(error.ne.0)then
   write(sparse%unlog,'(2(a,i0))')' The following error for phase ',phase,' was detected: ',error
   stop
  endif
 end subroutine

end subroutine
#else
subroutine solve_csr(sparse,x,y)
 !sparse*x=y
 class(csr),intent(in)::sparse
 real(kind=real8),intent(out)::x(:)
 real(kind=real8),intent(in)::y(:)

 write(sparse%unlog,'(a)')' Warning: Pardiso is not enabled!'
 x=y

end subroutine
#endif

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

subroutine printtobin_csr(sparse,un,typesp,namefile)
 class(csr),intent(in)::sparse
 integer(kind=int4),intent(in)::un
 character(len=*),intent(in),optional::typesp
 character(len=*),intent(in),optional::namefile

 integer(kind=int4)::i,unin
 integer(kind=intnel)::j
 character(len=10)::cdummy
 
 write(cdummy,'(i0)')un
 if(present(typesp).and..not.present(namefile))then
  open(newunit=unin,file='subpcg.'//adjustl(typesp(:len_trim(typesp)))//adjustl(cdummy(:len_trim(cdummy))),action='write',access='stream',buffered='yes')
 elseif(.not.present(typesp).and.present(namefile))then
  open(newunit=unin,file=adjustl(namefile(:len_trim(namefile)))//'.'//adjustl(cdummy(:len_trim(cdummy))),action='write',access='stream',buffered='yes')
 elseif(present(typesp).and.present(namefile))then
  open(newunit=unin,file=adjustl(typesp(:len_trim(typesp)))//'.'//&
                         adjustl(namefile(:len_trim(namefile)))//'.'//&
                         adjustl(cdummy(:len_trim(cdummy))),action='write',access='stream',buffered='yes')
 endif

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
 sparse%lsquare=.true.
 sparse%sizeia=0
 if(allocated(sparse%ia))deallocate(sparse%ia)
 if(allocated(sparse%ja))deallocate(sparse%ja)
 if(allocated(sparse%a))deallocate(sparse%a)

end subroutine

end module
