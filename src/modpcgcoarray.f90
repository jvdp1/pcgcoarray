module modpcgcoarray
 !$ use omp_lib
 !$ use mkl_service
 use modkind
 use modsparse
 implicit none
 private
 public::pcgrowcoarray

contains

!PUBLIC
subroutine pcgrowcoarray(crs,x,crhs,crsprecond,startrow,endrow,unlog)
 type(crssparse),intent(inout)::crs
 type(crssparse),intent(inout)::crsprecond
 integer(kind=int4),intent(in)::startrow[*],endrow[*]
 integer(kind=int4),intent(in)::unlog
 real(kind=real8),allocatable,intent(inout)::x(:)[:]
 character(len=*),intent(in)::crhs

 integer(kind=int4)::thisimage,unconv
 integer(kind=int4)::i,j,k,neq,nrow,iter
 integer(kind=int4)::startrowk,maxit=10000
 real(kind=real8)::tol=1.e-12
 real(kind=real8)::oldtau,conv,thr,beta
 real(kind=real8),allocatable::b_norm[:],resvec1[:],alpha[:],tau[:]
 real(kind=real8),allocatable::rhs(:)!,precond(:)
 real(kind=real8),allocatable::z(:),r(:),w(:)
 real(kind=real8),allocatable::p(:)[:]
 real(kind=real8),allocatable::T(:,:)
 !$ real(kind=real8)::t1,t2,val

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(unlog,'(/a,i0,a)')' Start image ',thisimage,' of the PCG Coarray subroutine...'

 write(unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel
 !$ write(unlog,'(" Number of threads for MKL   : ",i0)')mkl_get_max_threads() 

 sync all

 !read rhs
 allocate(rhs(crs%getdim(1)))
 call readrhs(rhs,startrow,endrow,unlog)

 write(unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage
 neq=crs%getdim(2)
 nrow=crs%getdim(1)

 !check the preconditioner
 write(unlog,'(/a,i0)')' Stats of the preconditioner for image ',thisimage
 call crsprecond%printstats()

 !invert preconditioner
 !do i=1,nrow
 ! if(precond(i).ne.0)precond(i)=1.d0/precond(i)
 !enddo
 
 !Initialistion PCG arrays
 allocate(z(nrow),w(nrow),r(nrow))
 allocate(p(neq)[*])
 r=0.d0;p=0.d0;z=0.d0;w=0.d0

 oldtau=1.d0

 !initiatilisation of r
 allocate(b_norm[*],resvec1[*],alpha[*],tau[*])
 if(.not.allocated(x))then
  allocate(x(neq)[*])
  x=0.d0
 endif
 if(sum(x).eq.0.d0)then
  r=rhs
 else
  print*,'not yet implemented'
  error stop
 endif

 b_norm=norm(rhs,1,nrow)
 resvec1=norm(r,1,nrow)

 sync all

 !update on all images
 !1. update on image 1
 if(thisimage.eq.1)then
  do i=2,num_images()
   !receives updates from other for its own image
   b_norm=b_norm+b_norm[i]
   resvec1=resvec1+resvec1[i]
  enddo
 endif
 sync all
 !2. update on the other image
 if(thisimage.ne.1)then
  b_norm=b_norm[1]
  resvec1=resvec1[1]
 endif
 sync all  !not sure if it is really needed

 conv=resvec1/b_norm

 write(unlog,'(/a,i0)')' Start iterating for image ',thisimage

 write(unlog,'(/a,e15.5)')' Norm of RHS: ',b_norm
 if(thisimage.eq.1)then
  open(newunit=unconv,file='convergence.dat',status='replace',action='write')
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')1,conv,0.
 endif

 if(thisimage.eq.num_images())then
  allocate(T(maxit,maxit));T=0.d0
 endif

 alpha=1.d0
 iter=2
 thr=tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.maxit)
  !z=Mi*r
  !do i=1,nrow
  ! z(i)=precond(i)*r(i)
  !enddo
  !M*z=r
  call crsprecond%solve(z,r)

  !tau=z*r
  tau=0.d0
  do i=1,nrow
   tau=tau+z(i)*r(i)
  enddo

  !update tau
  sync all
  if(thisimage.eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
    tau=tau+tau[i]
   enddo
  endif
  sync all
  !2. update on the other image
  if(thisimage.ne.1)tau=tau[1]

  beta=tau/oldtau
  oldtau=tau

  !p=z+beta*p
  do i=startrow,endrow
   p(i)=z(i-startrowk)+beta*p(i)
  enddo

  sync all
  do i=1,num_images()
   if(i.ne.thisimage)then
    !receives updates from other for its own image
    j=startrow[i]
    k=endrow[i]
    p(j:k)=p(j:k)[i]
   endif
  enddo
  sync all  !not sure if it is really needed

  !w=LHS*p
  call crs%multbyv(1._real8,'n',p,0._real8,w)
 
  !alpha=p*w
  if(thisimage.eq.num_images().and.iter.gt.2)call addalphabetatot(T,iter,alpha,beta)

  alpha=0.d0
  do i=1,nrow
   alpha=alpha+p(startrowk+i)*w(i)
  enddo

  !update alpha
  sync all
  if(thisimage.eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
    alpha=alpha+alpha[i]
   enddo
   alpha=tau/alpha
  endif
  sync all
  !2. update on the other image
  if(thisimage.ne.1)alpha=alpha[1]

  do i=startrow,endrow
   x(i)=x(i)+alpha*p(i)
  enddo

  do i=1,nrow
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,nrow)

  !update resvec1
  sync all
  if(thisimage.eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
    resvec1=resvec1+resvec1[i]
   enddo
  endif
  sync all
  !2. update on the other image
  if(thisimage.ne.1)resvec1=resvec1[1]

  conv=resvec1/b_norm

 if(thisimage.eq.1)then
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')iter,conv,resvec1
 endif

 iter=iter+1

 enddo
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(iter-1)
 !$ endif

 if(thisimage.eq.1)close(unconv)

 if(thisimage.eq.num_images())then
  !Estimation of eigenvalues and condition number
  if(iter>3)call eigenvalandcondnumber(T,iter,unlog)
  deallocate(T)
 endif

 sync all

 if(thisimage.eq.1)then
  do i=2,num_images()
   !receives updates from other for its own image
   j=startrow[i]
   k=endrow[i]
   x(j:k)=x(j:k)+x(j:k)[i]
  enddo
 endif
 sync all     !needed to wait for the update   ==>  probably better to use sync images

 write(unlog,'(/a,i0,a)')' End image ',thisimage,' of the PCG Coarray subroutine...'
 !$ write(unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2

end subroutine


!PRIVATE
!READ RHS
subroutine readrhs(rhs,startrow,endrow,unlog)
 !stupid to read a vector like that, but it works
 integer(kind=int4),intent(in)::startrow,endrow,unlog
 real(kind=real8),intent(inout)::rhs(:)

 integer(kind=int4)::i,j,io,un
 real(kind=real8)::val

 write(unlog,'(/a)')' Starts to read the rhs'
 rhs=0.d0
 open(newunit=un,file='rhs.bin',access='stream',action='read',status='old')!,buffered='yes')
 read(un)i
 i=0
 j=0
 do
  read(un,iostat=io)val
  if(io.ne.0)exit
  i=i+1
  if(i.gt.endrow)exit
  if(i.ge.startrow)then !.and.i.le.endrow)then
   j=j+1
   rhs(j)=val
  endif
 enddo
 close(un)
 write(unlog,'(a)')' End of reading the rhs'

end subroutine

!NORM
function norm(vector,starteq,endeq)
 real(kind=real8),intent(in)::vector(:)
 integer(kind=int4),intent(in)::starteq,endeq
 real(kind=real8)::norm
 
 integer(kind=int4)::i

 norm=0.d0
 do i=starteq,endeq
  norm=norm+vector(i)**2
 enddo

end function

!Condition number
subroutine addalphabetatot(T,iter,previousalpha,beta)
 integer(kind=int4),intent(in)::iter
 real(kind=real8),intent(in)::previousalpha,beta
 real(kind=real8),intent(inout)::T(:,:)
 
 integer(kind=int4)::previousiter
 real(kind=real8)::b

 b=sqrt(beta)
 previousiter=iter-1

 T(previousiter,previousiter)=T(previousiter,previousiter)+1.d0/previousalpha
 T(previousiter,iter)=T(previousiter,iter)+b/previousalpha
 T(iter,previousiter)=T(iter,previousiter)+b/previousalpha
 T(iter,iter)=T(iter,iter)+beta/previousalpha

end subroutine

subroutine eigenvalandcondnumber(T,iter,unlog)
 integer(kind=int4),intent(in)::iter
 integer(kind=int4),intent(in),optional::unlog
 real(kind=real8),intent(inout)::T(:,:)

 integer(kind=int4)::un
 real(kind=real8),allocatable::eigenvalue(:),ttmp(:,:)

 un=6
 if(present(unlog))un=unlog

 allocate(eigenvalue(iter-3))

 !call eigensyevd(T(2:iter-2,2:iter-2),iter-3,eigenvalue)

 allocate(ttmp,source=T(2:iter-2,2:iter-2))
 call eigensyevd(ttmp,iter-3,eigenvalue)
 deallocate(ttmp)

 write(un,'(/"  Mininum | maximum eigenvalues        : ",g0.5,a,g0.5)')minval(eigenvalue),' | ',maxval(eigenvalue)
 write(un,'( "   Effective spectral condition number : ",g0.5)')maxval(eigenvalue)/minval(eigenvalue)
 deallocate(eigenvalue)

end subroutine

subroutine eigensyevd(mat,n,w,leigvectors)
 integer(kind=int4),intent(in)::n
 real(kind=real8),intent(inout)::mat(:,:),w(:)
 logical,intent(in),optional::leigvectors

 integer(kind=int4)::info,i,j,lwork,liwork
 integer(kind=int4),allocatable::iwork(:)
 real(kind=real8),allocatable::work(:)
 character(len=1)::jobz

 jobz='N'
 if(present(leigvectors))then
  if(leigvectors)jobz='V'
 endif

 lwork=2*n*n+6*n+1
 liwork=5*n+3
 allocate(work(lwork),iwork(liwork))

 lwork=-1
 call dsyevd(jobz,'L',n,mat,n,w,work,lwork,iwork,liwork,info)
 if (info.eq.0) then
  lwork=int(work(1))
  liwork=iwork(1)
  deallocate(work,iwork)
  allocate(work(lwork),iwork(liwork))
 else
  print*,'ERROR info : ',info
  return
 endif

 call dsyevd(jobz,'L',n,mat,n,w,work,lwork,iwork,liwork,info)
 if (info.ne.0) then
  print*,'ERROR info-end: ',info
 endif

end subroutine

end module