!ifort -O3 -heap-arrays -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -parallel -lpthread -coarray=distributed -coarray-num-images=3 pcgcoarray.f90 -o pcgcoarray
!include 'modfiles.f90'

program  pcgcorray
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::io,un,i,j,k,l,m,n,neq,ncol,iter
 integer(kind=int4)::startcolk,maxit=2000
 integer(kind=int4)::startrow[*],endrow[*],startcol[*],endcol[*]
 integer(kind=intnel)::nel
 character(len=80)::host,cdummy
 real(kind=real8)::tol=1.e-12
 real(kind=real8)::oldtau,conv,thr,beta
 real(kind=real8)::b_norm[*],resvec1[*],alpha[*],tau[*]
 real(kind=real8),allocatable::rhs(:),precond(:)
 real(kind=real8),allocatable::p(:),z(:)
 real(kind=real8),allocatable::r(:)[:],w(:)[:],x(:)[:]
 real(kind=real8)::val
 type(csr)::sparse

 call get_environment_variable("HOSTNAME",value=host)
 write(*,*)"Hello from image ",this_image()," out of ",num_images()," total images on host",trim(host)

 sync all

 if(this_image().eq.1)then
  open(newunit=un,file='param.pcgcoarray',status='old',action='read')
  n=0
  do
   read(un,*,iostat=io)i,j,k,l,m
   if(io.ne.0)exit
   n=n+1
   if(n.gt.num_images())then
    write(*,'(a)')' The parameter file is not correct!'
    stop
   endif
   startrow[i]=j
   endrow[i]=k
   startcol[i]=l
   endcol[i]=m
  enddo
  close(un)
  if(n.ne.num_images())then
   write(*,'(a)')' The parameter file is not correct!'
   stop
  endif
 endif

 sync all

 !Reads the matrix
 write(*,'(/a,i0)')' Reads the matrix from image ',this_image()
 write(cdummy,'(i0)')this_image()

 open(newunit=un,file='subpcg.'//adjustl(cdummy(:len_trim(cdummy))),status='old',action='read',access='stream')
 read(un)i,j,nel
 call sparse%alloc(nel,i,j)
 
 read(un)sparse%ia
 print*,'read ia'
 read(un)sparse%ja
 print*,'read ja'
 read(un)sparse%a
 print*,'read a'
 close(un)


 write(*,'(a,i0/)')' End of reading the matrix from image ',this_image()

 !call sparse%printfile(600+this_image())

 !sync all

 !read rhs
 write(*,*)"Image ",this_image()," starts to read the rhs"
 allocate(rhs(sparse%m))
 rhs=0.d0
 open(newunit=un,file='rhs.bin',access='stream',action='read',status='old',buffered='yes')
 read(un)i
 i=0
 j=0
 do
  read(un,iostat=io)val
  if(io.ne.0)exit
  i=i+1
  if(i.gt.endcol)exit
  if(i.ge.startcol)then !.and.i.le.endcol)then
   j=j+1
   rhs(j)=val
  endif
 enddo
 close(un)
 !write(*,*)"Hello from image ",this_image(),'aaa',rhs

 !sync all


 write(*,*)' Preparation for the PCG...'
 neq=sparse%n
 ncol=sparse%m
 !create preconditioner
 precond=sparse%diag(startcol)
 do i=1,ncol
  if(precond(i).ne.0)precond(i)=1.d0/precond(i)
 enddo
 
 !Initialistion PCG arrays
 allocate(x(neq)[*])
 x=0.d0

 allocate(p(ncol),z(ncol))
 allocate(r(neq)[*],w(neq)[*])
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

 sync all

! !update of r
! do i=1,num_images()
!  !receives updates from other for its own image
!  if(i.ne.this_image())then
!   j=startcol[i]
!   k=endcol[i]
!   r(j:k)=r(j:k)+r(j:k)[i]
!   b_norm=b_norm+b_norm[i]
!  endif
! enddo

 !update on all images
 !1. update on image 1
 if(this_image().eq.1)then
  do i=2,num_images()
   !receives updates from other for its own image
   j=startcol[i]
   k=endcol[i]
   r(j:k)=r(j:k)+r(j:k)[i]
   b_norm=b_norm+b_norm[i]
   resvec1=resvec1+resvec1[i]
  enddo
 endif
 sync all
 !2. update on the other image
 if(this_image().ne.1)then
  r(:)=r(:)[1]
  b_norm=b_norm[1]
  resvec1=resvec1[1]
 endif
 sync all  !not sure if it is needed

 conv=resvec1/b_norm

 if(this_image().eq.1)then
  write(*,'(a,e15.5)')' Norm of RHS: ',b_norm
  write(*,'(" Iteration ",i6," Convergence = ",e12.5)')1,conv
 endif

 alpha=1.d0
 iter=2
 thr=tol*b_norm

 !Start iteration
 startcolk=startcol-1
 do while(resvec1.gt.thr.and.iter.le.maxit)
  !z=M*r
  do i=1,ncol
   z(i)=precond(i)*r(startcolk+i)
  enddo

  !tau=z*r
  tau=0.d0
  do i=1,ncol
   tau=tau+z(i)*r(startcolk+i)
  enddo

  !update tau
  sync all
  if(this_image().eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
    tau=tau+tau[i]
   enddo
  endif
  sync all
  !2. update on the other image
  if(this_image().ne.1)tau=tau[1]

  beta=tau/oldtau
  oldtau=tau

  !p=z+beta*p
  do i=1,ncol
   p(i)=z(i)+beta*p(i)
  enddo

  !w=LHS*p
  call multgenv(sparse,p,w)
 
  !update w
  sync all
  if(this_image().eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
   w=w+w(:)[i]
   enddo
  endif
  sync all
  !2. update on the other image
  if(this_image().ne.1)w(:)=w(:)[1]

  !alpha=p*w
  alpha=0.d0
  do i=1,ncol
   alpha=alpha+p(i)*w(startcolk+i)
  enddo

  !update alpha
  sync all
  if(this_image().eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
    alpha=alpha+alpha[i]
   enddo
   alpha=tau/alpha
  endif
  sync all
  !2. update on the other image
  if(this_image().ne.1)alpha=alpha[1]

  do i=1,ncol
   x(startcolk+i)=x(startcolk+i)+alpha*p(i)
  enddo

!  do i=startcolk+1,startcolk+ncol
!   r(i)=r(i)-alpha*w(i)
!  enddo
!  resvec1=norm(r,startcolk+1,startcolk+ncol)
!  !update
!  sync all
!  do i=1,num_images()
!   if(i.ne.this_image())then
!    j=startcol[i]
!    k=endcol[i]
!    r(j:k)=r(j:k)[i]
!    resvec1=resvec1+resvec1[i]
!   endif
!  enddo

  do i=1,neq
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,neq)

  conv=resvec1/b_norm

!  sync all

  if(this_image().eq.1)then
   write(*,'(" Iteration ",i6," Convergence = ",e12.5,x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo

 sync all
 if(this_image().eq.1)then
  do i=2,num_images()
   x=x+x(:)[i]
  enddo
  call print_ascii(x,1,neq)
 endif


! sync all 
! if(this_image().eq.1)then
!   write(*,*)i,'aaa',size(r(:)),'ccc',r(:)
!   write(*,*)i,'aaa',b_norm
!   write(*,*)i,'aaa',x
!  do i=1,num_images()
!   write(*,*)i,'aaa',size(r(:)[i]),'ccc',r(:)[i]
!   write(*,*)i,'aaa',b_norm
!  enddo
! endif

 sync all
 write(*,*)"End for image ",this_image()," out of ",num_images()," total images!"

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

 open(newunit=un,file='solutions.pcgcoarray')
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
