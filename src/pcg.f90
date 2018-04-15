!ifort -O3 -heap-arrays -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -parallel -lpthread  pcg.f90 -o pcg

program  pcg
 !$ use omp_lib
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::io,un,i,j,k,l,m,n,neq,ncol,iter
 integer(kind=int4)::startcolk,maxit=4000
 integer(kind=int4)::startrow,endrow,startcol,endcol
 integer(kind=intnel)::nel
 character(len=80)::host,cdummy
 real(kind=real8)::tol=1.e-12
 real(kind=real8)::oldtau,conv,thr,beta
 real(kind=real8)::b_norm,resvec1,alpha,tau
 real(kind=real8),allocatable::rhs(:),precond(:)
 real(kind=real8),allocatable::x(:),p(:),z(:)
 real(kind=real8),allocatable::r(:),w(:)
 real(kind=real8)::t1,val
 type(csr)::sparse


 startcol=1
 !Reads the matrix

 open(newunit=un,file='matrixija.bin',status='old',action='read',access='stream')
 read(un)i,nel
 call sparse%alloc(nel,i-1)
 
 read(un)sparse%ia
 print*,'read ia'
 read(un)sparse%ja
 print*,'read ja'
 read(un)sparse%a
 print*,'read a'
 close(un)
 write(*,'(a,i0/)')' End of reading the matrix '


 !read rhs
 write(*,*)"read the rhs"
 allocate(rhs(sparse%m))
 rhs=0.d0
 open(newunit=un,file='rhs.bin',access='stream',action='read',status='old',buffered='yes')
 read(un)i
 read(un)rhs

 write(*,*)' Preparation for the PCG...'
 neq=sparse%n
 ncol=sparse%m
 endcol=sparse%m
 !create preconditioner
 precond=sparse%diagcol(startcol)
 do i=1,ncol
  if(precond(i).ne.0)precond(i)=1.d0/precond(i)
 enddo
 
 !Initialistion PCG arrays
 allocate(x(neq))
 x=0.d0

 allocate(p(ncol),z(ncol))
 allocate(r(neq),w(neq))
 r=0.d0;p=0.d0;z=0.d0;w=0.d0

 oldtau=1.d0

 !initiatilisation of r
 if(sum(x).eq.0.d0)then
  r(startcol:endcol)=rhs
 else
  print*,'not yet implemented'
  stop
 endif

 b_norm=norm(rhs,1,sparse%m)
 resvec1=norm(r,startcol,endcol)

 conv=resvec1/b_norm

  write(*,'(a,e15.5)')' Norm of RHS: ',b_norm
  write(*,'(" Iteration ",i6," Convergence = ",e12.5)')1,conv

 alpha=1.d0
 iter=2
 thr=tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 do while(resvec1.gt.thr.and.iter.le.maxit)
  startcolk=startcol-1
  !z=M*r
  do i=1,ncol
   z(i)=precond(i)*r(startcolk+i)
  enddo

  !tau=z*r
  tau=0.d0
  do i=1,ncol
   tau=tau+z(i)*r(startcolk+i)
  enddo

!print*,'aaa',tau
  beta=tau/oldtau
  oldtau=tau

  !p=z+beta*p
  do i=1,ncol
   p(i)=z(i)+beta*p(i)
  enddo
!print*,'bbb',tau,beta,p
!stop

  !w=LHS*p
  call multgenv(sparse,p,w)
!print*,'bbb',tau,beta,w
!stop

  !alpha=p*w
  alpha=0.d0
  do i=1,ncol
   alpha=alpha+p(i)*w(startcolk+i)
  enddo

  alpha=tau/alpha
!print*,'ddd',alpha

  do i=1,ncol
   x(startcolk+i)=x(startcolk+i)+alpha*p(i)
  enddo

  do i=startcolk+1,startcolk+ncol
   r(i)=r(i)-alpha*w(i)
  enddo
 
  resvec1=norm(r,startcolk+1,startcolk+ncol)

  conv=resvec1/b_norm

  write(*,'(" Iteration ",i6," Convergence = ",e12.5,x,e12.5)')iter,conv,resvec1
  iter=iter+1
 enddo

 !$ val=omp_get_wtime()-t1
 !$  write(*,'("  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(*,'("  Approximate Wall clock time per iteration (seconds): ",f12.2)')val/(iter-1)
 
 call print_ascii(x,1,neq)

!   write(*,*)i,'aaa',size(r(:)),'ccc',r(:)
!   write(*,*)i,'aaa',b_norm
!   write(*,*)i,'aaa',x


contains

subroutine multgenv(A,x,y)
 !Computes y=0.d0*y+A*x
 type(csr),intent(in)::A
 real(kind=real8),intent(in)::x(:)
 real(kind=real8),intent(out)::y(:)

 character(len=1)::matdescra(6)

 matdescra(1)='G'
 matdescra(4)='F'
 
 call mkl_dcsrmv('N',A%n,A%m,1.d0,matdescra,A%a,A%ja,A%ia(1:A%n),A%ia(2:A%n+1),x,0.d0,y)
 
end subroutine

subroutine print_ascii(x,startpos,endpos)
 integer(kind=intc),intent(in)::startpos,endpos
 real(kind=real8),intent(in)::x(:)

 integer(kind=intc)::un

 open(newunit=un,file='solutions.pcg')
 do i=startpos,endpos
  write(un,*)x(i)
 enddo
 close(un)
end subroutine

function norm(vector,starteq,endeq)
 real(kind=realc),intent(in)::vector(:)
 integer(kind=intc),intent(in)::starteq,endeq
 real(kind=realc)::norm
 
 integer(kind=intc)::i

 norm=0.d0
 do i=starteq,endeq
  norm=norm+vector(i)**2
 enddo

end function

end program
