module modcoarraysolver
 use iso_fortran_env
 !$ use omp_lib
 use modsparse
 use modprecond,only:gen_precond
 use modcoeff,only:gen_coeff
 use modcoarrayop
 implicit none
 private
 public::cg,chebyshev,pcg,fpcg

 !ABSTRACTS
 type,abstract::solver
  private
  integer(kind=int32)::neq=-9
  integer(kind=int32)::unlog=6
  integer(kind=int32)::maxit=10000
  integer(kind=int32)::niter=-1
  real(kind=real64)::tol=1.e-6
  contains
  private
  procedure,public::getiter
  procedure,public::setmaxiterations
  procedure,public::setthreshold
  procedure,public::setoutput
 end type

 type,extends(solver),abstract::cgsolver
  private
  real(kind=real64)::smalleigenval=0._real64
  real(kind=real64)::largeeigenval=0._real64
  contains
  private
  procedure,public::getsmalleigenvalue
  procedure,public::getlargeeigenvalue
 end type


 !---SOLVERS
 !*CHEBYSHEV
 type,extends(solver)::chebyshev
  private
  real(kind=real64)::smalleigenval,largeeigenval
  contains
  private
  procedure,public::solve=>chebyshevrowcoarray
  final::destroy_chebyshev
 end type
 
 interface chebyshev
  module procedure constructor_chebyshev
 end interface

 
 !*CG
 type,extends(cgsolver)::cg
  private
  contains
  private
  procedure,public::solve=>cgrowcoarray
  final::destroy_cg
 end type
 
 interface cg
  module procedure constructor_cg
 end interface


 !*PCG
 type,extends(cgsolver)::pcg
  private
  contains
  private
  procedure,public::solve=>pcgrowcoarray
  final::destroy_pcg
 end type
 
 interface pcg
  module procedure constructor_pcg
 end interface


 !*FLEXIBLE PCG
 type,extends(cgsolver)::fpcg
  private
  contains
  private
  procedure,public::solve=>fpcgrowcoarray
  final::destroy_fpcg
 end type
 
 interface fpcg
  module procedure constructor_fpcg
 end interface

contains

!PUBLIC

!--ABSTRACT
!**GET NUMBER ITERATIONS
function getiter(this) result(niter)
 class(solver),intent(in)::this
 integer(kind=int32)::niter

 niter=this%niter

end function

!**SET MAXIMUM ITERATIONS
subroutine setmaxiterations(this,maxit)
 class(solver),intent(inout)::this
 integer(kind=int32),intent(in)::maxit

 this%maxit=maxit

end subroutine

!**SET OUTPUT
subroutine setoutput(this,un)
 class(solver),intent(inout)::this
 integer(kind=int32),intent(in)::un

 this%unlog=un

end subroutine

!**SET THRESHOLD
subroutine setthreshold(this,tol)
 class(solver),intent(inout)::this
 real(kind=real64),intent(in)::tol

 this%tol=tol

end subroutine


!--CHEBYSHEV
!**CONSTRUCTOR
function constructor_chebyshev(neq,small,large) result(this)
 type(chebyshev)::this
 integer(kind=int32),intent(in)::neq
 real(kind=real64),intent(in)::small,large

 this%neq=neq
 
 this%smalleigenval=small
 this%largeeigenval=large

end function

!**SET EIGENVALUES
subroutine seteigenvalues(this,small,large)
 class(chebyshev),intent(inout)::this
 real(kind=real64),intent(in)::small,large

 this%smalleigenval=small
 this%largeeigenval=large

end subroutine

!**SOLVER
subroutine chebyshevrowcoarray(this,crs,x,crhs,precond,startrow,endrow)
 class(chebyshev),intent(inout)::this
 class(*),intent(inout)::crs
 class(*),intent(inout)::precond
 integer(kind=int32),intent(in)::startrow[*],endrow[*]
 real(kind=real64),intent(inout)::x(:)[*]
 character(len=*),intent(in)::crhs

 integer(kind=int32)::thisimage,unconv
 integer(kind=int32)::i,j,k,nrow,iter
 integer(kind=int32)::startrowk
 real(kind=real64)::conv,thr,beta
 real(kind=real64)::cvalue,dvalue
 real(kind=real64),allocatable::b_norm[:],resvec1[:],alpha[:]
 real(kind=real64),allocatable::rhs(:)
 real(kind=real64),allocatable::z(:),r(:),w(:)
 real(kind=real64),allocatable::p(:)[:]
 !$ real(kind=real64)::t1,t2,val

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(this%unlog,'(//a,i0,a)')' **Start image ',thisimage,' of the CHEBYSHEV Coarray subroutine...'

 write(this%unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(this%unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel

 !check provided eigenvalues
 if(this%smalleigenval.eq.0._real64.or.this%largeeigenval.eq.0._real64)then
  write(this%unlog,'(a)')' ERROR: Wrong eigenvalues for the CHEBYSHEV solver!'
  error stop
 endif

 sync all

 nrow=endrow-startrow+1

 !read rhs
 allocate(rhs(nrow))
 call readrhs(rhs,crhs,startrow,endrow,this%unlog)

 write(this%unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage

 !Initialistion PCG arrays
 allocate(z(nrow),w(nrow),r(nrow))
 allocate(p(this%neq)[*])
 r=0._real64;p=0._real64;z=0._real64;w=0._real64

 !initiatilisation of r
 allocate(b_norm[*],resvec1[*],alpha[*])

 if(sum(x).eq.0._real64)then
  r=rhs
 else
  print*,'not yet implemented'
  error stop
 endif

 b_norm=norm(rhs,1,nrow)
 resvec1=norm(r,1,nrow)

 sync all

 !update on all images
! !1. update on image 1
! if(thisimage.eq.1)then
!  do i=2,num_images()
!   !receives updates from other for its own image
!   b_norm=b_norm+b_norm[i]
!   resvec1=resvec1+resvec1[i]
!  enddo
! endif
! sync all
! !2. update on the other image
! if(thisimage.ne.1)then
!  b_norm=b_norm[1]
!  resvec1=resvec1[1]
! endif
! sync all  !not sure if it is really needed
 call co_sum(b_norm)
 call co_sum(resvec1)

 conv=resvec1/b_norm

 write(this%unlog,'(/a,i0)')' Start iterating for image ',thisimage

 write(this%unlog,'(/a,e15.5)')' Norm of RHS: ',b_norm
 if(thisimage.eq.1)then
  open(newunit=unconv,file='convergence.dat',status='replace',action='write')
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')1,conv,0.
 endif

 dvalue=(this%largeeigenval+this%smalleigenval)/2._real64
 cvalue=(this%largeeigenval-this%smalleigenval)/2._real64

 alpha=0._real64
 iter=2
 thr=this%tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.this%maxit)
  !M*z=r
  call solveprecond(precond,z,r,this%unlog)

  beta=(cvalue*alpha/2._real64)**2
  alpha=1._real64/(dvalue-beta)

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

  do i=startrow,endrow
   x(i)=x(i)+alpha*p(i)
  enddo

  !w=LHS*p
  call multbyv(crs,p,w,startrow,endrow,this%unlog)

  do i=1,nrow
   r(i)=r(i)-alpha*w(i)
  enddo

  resvec1=norm(r,1,nrow)

  !update resvec1
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    resvec1=resvec1+resvec1[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)resvec1=resvec1[1]
  call co_sum(resvec1)

  conv=resvec1/b_norm

  if(thisimage.eq.1)then
   write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 this%niter=iter-1
 write(this%unlog,'(/"  Number of iterations                               : ",i0)')this%niter
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(this%unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(this%unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(this%niter)
 !$ endif

 if(thisimage.eq.1)close(unconv)

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

 write(this%unlog,'(/a,i0,a)')' **End image ',thisimage,' of the PCG Coarray subroutine...'
 !$ write(this%unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2

end subroutine

!*FINAL
subroutine destroy_chebyshev(this)
 type(chebyshev),intent(inout)::this

 call destroy_solver(this)

 this%smalleigenval=0._real64
 this%largeeigenval=0._real64
 
end subroutine


!--CG SOLVER ABSTRACT
!**GET EIGENVALUES
function getlargeeigenvalue(this) result(val)
 class(cgsolver),intent(in)::this
 real(kind=real64)::val
 
 val=this%largeeigenval

end function

function getsmalleigenvalue(this) result(val)
 class(cgsolver),intent(in)::this
 real(kind=real64)::val
 
 val=this%smalleigenval

end function


!--CG
!**CONSTRUCTOR
function constructor_cg(neq) result(this)
 type(cg)::this
 integer(kind=int32),intent(in)::neq

 this%neq=neq

end function

!**SOLVER
subroutine cgrowcoarray(this,crs,x,crhs,startrow,endrow)
 class(cg),intent(inout)::this
 class(*),intent(inout)::crs
 integer(kind=int32),intent(in)::startrow[*],endrow[*]
 real(kind=real64),intent(inout)::x(:)[*]
 character(len=*),intent(in)::crhs

 integer(kind=int32)::thisimage,unconv
 integer(kind=int32)::i,j,k,nrow,iter
 integer(kind=int32)::startrowk
 real(kind=real64)::oldtau,conv,thr,beta
 real(kind=real64),allocatable::b_norm[:],resvec1[:],alpha[:],tau[:]
 real(kind=real64),allocatable::rhs(:)
 real(kind=real64),allocatable::r(:),w(:)
 real(kind=real64),allocatable::p(:)[:]
 real(kind=real64),allocatable::T(:,:)
 !$ real(kind=real64)::t1,t2,val

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(this%unlog,'(//a,i0,a)')' **Start image ',thisimage,' of the PCG Coarray subroutine...'

 write(this%unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(this%unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel

 sync all

 nrow=endrow-startrow+1

 !read rhs
 allocate(rhs(nrow))
 call readrhs(rhs,crhs,startrow,endrow,this%unlog)

 write(this%unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage

 !Initialistion PCG arrays
 allocate(w(nrow),r(nrow))
 allocate(p(this%neq)[*])
 r=0._real64;p=0._real64;w=0._real64

 oldtau=1._real64

 !initiatilisation of r
 allocate(b_norm[*],resvec1[*],alpha[*],tau[*])

 if(sum(x).eq.0._real64)then
  r=rhs
 else
  print*,'not yet implemented'
  error stop
 endif

 b_norm=norm(rhs,1,nrow)
 resvec1=norm(r,1,nrow)

 sync all

 !update on all images
! !1. update on image 1
! if(thisimage.eq.1)then
!  do i=2,num_images()
!   !receives updates from other for its own image
!   b_norm=b_norm+b_norm[i]
!   resvec1=resvec1+resvec1[i]
!  enddo
! endif
! sync all
! !2. update on the other image
! if(thisimage.ne.1)then
!  b_norm=b_norm[1]
!  resvec1=resvec1[1]
! endif
 !sync all  !not sure if it is really needed
 call co_sum(b_norm)
 call co_sum(resvec1)

 conv=resvec1/b_norm

 write(this%unlog,'(/a,i0)')' Start iterating for image ',thisimage

 write(this%unlog,'(/a,e15.5)')' Norm of RHS: ',b_norm
 if(thisimage.eq.1)then
  open(newunit=unconv,file='convergence.dat',status='replace',action='write')
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')1,conv,0.
 endif

 if(thisimage.eq.num_images())then
  allocate(T(this%maxit,this%maxit));T=0._real64
 endif

 alpha=1._real64
 iter=2
 thr=this%tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.this%maxit)

  !tau=r*r
  tau=0._real64
  do i=1,nrow
   tau=tau+r(i)**2
  enddo

  !update tau
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    tau=tau+tau[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)tau=tau[1]
  call co_sum(tau)

  beta=tau/oldtau
  oldtau=tau

  !p=r+beta*p
  do i=startrow,endrow
   p(i)=r(i-startrowk)+beta*p(i)
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
  call multbyv(crs,p,w,startrow,endrow,this%unlog)
 
  !alpha=p*w
  if(thisimage.eq.num_images().and.iter.gt.2)call addalphabetatot(T,iter,alpha,beta)

  alpha=0._real64
  do i=1,nrow
   alpha=alpha+p(startrowk+i)*w(i)
  enddo

  !update alpha
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    alpha=alpha+alpha[i]
!   enddo
!   alpha=tau/alpha
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)alpha=alpha[1]
  call co_sum(alpha)
  alpha=tau/alpha

  do i=startrow,endrow
   x(i)=x(i)+alpha*p(i)
  enddo

  do i=1,nrow
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,nrow)

  !update resvec1
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    resvec1=resvec1+resvec1[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)resvec1=resvec1[1]
  call co_sum(resvec1)

  conv=resvec1/b_norm

  if(thisimage.eq.1)then
   write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 this%niter=iter-1
 write(this%unlog,'(/"  Number of iterations                               : ",i0)')this%niter
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(this%unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(this%unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(this%niter)
 !$ endif

 if(thisimage.eq.1)close(unconv)

 if(thisimage.eq.num_images())then
  !Estimation of eigenvalues and condition number
  if(iter>3)call eigenvalandcondnumber(T,this%niter,this%smalleigenval,this%largeeigenval,this%unlog)
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
! if(thisimage.ne.1)then
!   j=startrow
!   k=endrow
!   x(j:k)[1]=x(j:k)[1]+x(j:k)
! endif

 sync all     !needed to wait for the update   ==>  probably better to use sync images

 write(this%unlog,'(/a,i0,a)')' **End image ',thisimage,' of the PCG Coarray subroutine...'
 !$ write(this%unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2

end subroutine

!*FINAL
subroutine destroy_cg(this)
 type(cg),intent(inout)::this

 call destroy_solver(this)

 this%smalleigenval=0._real64
 this%largeeigenval=0._real64

end subroutine


!--PCG
!**CONSTRUCTOR
function constructor_pcg(neq) result(this)
 type(pcg)::this
 integer(kind=int32),intent(in)::neq

 this%neq=neq

end function

!**SOLVER
subroutine pcgrowcoarray(this,crs,x,crhs,precond,startrow,endrow)
 class(pcg),intent(inout)::this
 class(*),intent(inout)::crs
 class(*),intent(inout),optional::precond
 integer(kind=int32),intent(in)::startrow[*],endrow[*]
 real(kind=real64),intent(inout)::x(:)[*]
 character(len=*),intent(in)::crhs

 integer(kind=int32)::thisimage,unconv
 integer(kind=int32)::i,j,k,nrow,iter
 integer(kind=int32)::startrowk
 real(kind=real64)::oldtau,conv,thr,beta
 real(kind=real64),allocatable::b_norm[:],resvec1[:],alpha[:],tau[:]
 real(kind=real64),allocatable::rhs(:)!,precond(:)
 real(kind=real64),allocatable::z(:),r(:),w(:)
 real(kind=real64),allocatable::p(:)[:]
 real(kind=real64),allocatable::T(:,:)
 !$ real(kind=real64)::t1,t2,val

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(this%unlog,'(//a,i0,a)')' **Start image ',thisimage,' of the PCG Coarray subroutine...'

 write(this%unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(this%unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel

 sync all

 nrow=endrow-startrow+1

 !read rhs
 allocate(rhs(nrow))
 call readrhs(rhs,crhs,startrow,endrow,this%unlog)

 write(this%unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage

 !Initialistion PCG arrays
 allocate(z(nrow),w(nrow),r(nrow))
 allocate(p(this%neq)[*])
 r=0._real64;p=0._real64;z=0._real64;w=0._real64

 oldtau=1._real64

 !initiatilisation of r
 allocate(b_norm[*],resvec1[*],alpha[*],tau[*])

 if(sum(x).eq.0._real64)then
  r=rhs
 else
  print*,'not yet implemented'
  error stop
 endif

 b_norm=norm(rhs,1,nrow)
 resvec1=norm(r,1,nrow)

 sync all

 !update on all images
! !1. update on image 1
! if(thisimage.eq.1)then
!  do i=2,num_images()
!   !receives updates from other for its own image
!   b_norm=b_norm+b_norm[i]
!   resvec1=resvec1+resvec1[i]
!  enddo
! endif
! sync all
! !2. update on the other image
! if(thisimage.ne.1)then
!  b_norm=b_norm[1]
!  resvec1=resvec1[1]
! endif
 call co_sum(b_norm)
 call co_sum(resvec1)

 conv=resvec1/b_norm

 write(this%unlog,'(/a,i0)')' Start iterating for image ',thisimage

 write(this%unlog,'(/a,e15.5)')' Norm of RHS: ',b_norm
 if(thisimage.eq.1)then
  open(newunit=unconv,file='convergence.dat',status='replace',action='write')
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')1,conv,0.
 endif

 if(thisimage.eq.num_images())then
  allocate(T(this%maxit,this%maxit));T=0._real64
 endif

 alpha=1._real64
 iter=2
 thr=this%tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.this%maxit)
  !M*z=r
  call solveprecond(precond,z,r,this%unlog)

  !tau=z*r
  tau=0._real64
  do i=1,nrow
   tau=tau+z(i)*r(i)
  enddo

  !update tau
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    tau=tau+tau[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)tau=tau[1]
  call co_sum(tau)

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
  call multbyv(crs,p,w,startrow,endrow,this%unlog)
 
  !alpha=p*w
  if(thisimage.eq.num_images().and.iter.gt.2)call addalphabetatot(T,iter,alpha,beta)

  alpha=0._real64
  do i=1,nrow
   alpha=alpha+p(startrowk+i)*w(i)
  enddo

  !update alpha
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    alpha=alpha+alpha[i]
!   enddo
!   alpha=tau/alpha
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)alpha=alpha[1]
  call co_sum(alpha)
  alpha=tau/alpha

  do i=startrow,endrow
   x(i)=x(i)+alpha*p(i)
  enddo

  do i=1,nrow
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,nrow)

  !update resvec1
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    resvec1=resvec1+resvec1[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)resvec1=resvec1[1]
  call co_sum(resvec1)

  conv=resvec1/b_norm

  if(thisimage.eq.1)then
   write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 this%niter=iter-1
 write(this%unlog,'(/"  Number of iterations                               : ",i0)')this%niter
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(this%unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(this%unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(this%niter)
 !$ endif

 if(thisimage.eq.1)close(unconv)

 if(thisimage.eq.num_images())then
  !Estimation of eigenvalues and condition number
  if(iter>3)call eigenvalandcondnumber(T,this%niter,this%smalleigenval,this%largeeigenval,this%unlog)
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
! if(thisimage.ne.1)then
!   j=startrow
!   k=endrow
!   x(j:k)[1]=x(j:k)[1]+x(j:k)
! endif

 sync all     !needed to wait for the update   ==>  probably better to use sync images

 write(this%unlog,'(/a,i0,a)')' **End image ',thisimage,' of the PCG Coarray subroutine...'
 !$ write(this%unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2

end subroutine

!*FINAL
subroutine destroy_pcg(this)
 type(pcg),intent(inout)::this

 call destroy_solver(this)

 this%smalleigenval=0._real64
 this%largeeigenval=0._real64

end subroutine

!--FLEXIBLE PCG
!**CONSTRUCTOR
function constructor_fpcg(neq) result(this)
 type(fpcg)::this
 integer(kind=int32),intent(in)::neq

 this%neq=neq
 
end function

!**SOLVER
subroutine fpcgrowcoarray(this,crs,x,crhs,precond,startrow,endrow)
 class(fpcg),intent(inout)::this
 class(*),intent(inout)::crs
 class(*),intent(inout),optional::precond
 integer(kind=int32),intent(in)::startrow[*],endrow[*]
 real(kind=real64),intent(inout)::x(:)[*]
 character(len=*),intent(in)::crhs

 integer(kind=int32)::thisimage,unconv
 integer(kind=int32)::i,j,k,nrow,iter
 integer(kind=int32)::startrowk
 real(kind=real64)::oldtau,conv,thr,beta
 real(kind=real64),allocatable::b_norm[:],resvec1[:],alpha[:],tau[:]
 real(kind=real64),allocatable::rhs(:)!,precond(:)
 real(kind=real64),allocatable::z(:),r(:),rp(:),w(:)
 real(kind=real64),allocatable::p(:)[:]
 real(kind=real64),allocatable::T(:,:)
 !$ real(kind=real64)::t1,t2,val

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(this%unlog,'(//a,i0,a)')' **Start image ',thisimage,' of the PCG Coarray subroutine...'

 write(this%unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(this%unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel

 sync all

 nrow=endrow-startrow+1

 !read rhs
 allocate(rhs(nrow))
 call readrhs(rhs,crhs,startrow,endrow,this%unlog)

 write(this%unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage

 !Initialistion PCG arrays
 allocate(z(nrow),w(nrow),r(nrow),rp(nrow))
 allocate(p(this%neq)[*])
 r=0._real64;rp=0._real64;p=0._real64;z=0._real64;w=0._real64

 oldtau=1._real64

 !initiatilisation of r
 allocate(b_norm[*],resvec1[*],alpha[*],tau[*])

 if(sum(x).eq.0._real64)then
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
! if(thisimage.eq.1)then
!  do i=2,num_images()
!   !receives updates from other for its own image
!   b_norm=b_norm+b_norm[i]
!   resvec1=resvec1+resvec1[i]
!  enddo
! endif
! sync all
! !2. update on the other image
! if(thisimage.ne.1)then
!  b_norm=b_norm[1]
!  resvec1=resvec1[1]
! endif
! sync all  !not sure if it is really needed
 call co_sum(b_norm)
 call co_sum(resvec1)

 conv=resvec1/b_norm

 write(this%unlog,'(/a,i0)')' Start iterating for image ',thisimage

 write(this%unlog,'(/a,e15.5)')' Norm of RHS: ',b_norm
 if(thisimage.eq.1)then
  open(newunit=unconv,file='convergence.dat',status='replace',action='write')
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')1,conv,0.
 endif

 if(thisimage.eq.num_images())then
  allocate(T(this%maxit,this%maxit));T=0._real64
 endif

 alpha=1._real64
 iter=2
 thr=this%tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.this%maxit)
  !M*z=r
  call solveprecond(precond,z,r,this%unlog)

  !tau=z*r
  tau=0._real64
  do i=1,nrow
   tau=tau+z(i)*(r(i)-rp(i))
  enddo

  !update tau
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    tau=tau+tau[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)tau=tau[1]
  call co_sum(tau)

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
  call multbyv(crs,p,w,startrow,endrow,this%unlog)
 
  !alpha=p*w
  if(thisimage.eq.num_images().and.iter.gt.2)call addalphabetatot(T,iter,alpha,beta)

  alpha=0._real64
  do i=1,nrow
   alpha=alpha+p(startrowk+i)*w(i)
  enddo

  !update alpha
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    alpha=alpha+alpha[i]
!   enddo
!   alpha=tau/alpha
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)alpha=alpha[1]
  call co_sum(alpha)
  alpha=tau/alpha

  do i=startrow,endrow
   x(i)=x(i)+alpha*p(i)
  enddo

  rp=r
  do i=1,nrow
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,nrow)

  !update resvec1
!  sync all
!  if(thisimage.eq.1)then
!   do i=2,num_images()
!    !receives updates from other for its own image
!    resvec1=resvec1+resvec1[i]
!   enddo
!  endif
!  sync all
!  !2. update on the other image
!  if(thisimage.ne.1)resvec1=resvec1[1]
  call co_sum(resvec1)

  conv=resvec1/b_norm

  if(thisimage.eq.1)then
   write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,1x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 this%niter=iter-1
 write(this%unlog,'(/"  Number of iterations                               : ",i0)')this%niter
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(this%unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(this%unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(this%niter)
 !$ endif

 if(thisimage.eq.1)close(unconv)

 if(thisimage.eq.num_images())then
  !Estimation of eigenvalues and condition number
  if(iter>3)call eigenvalandcondnumber(T,this%niter,this%smalleigenval,this%largeeigenval,this%unlog)
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
! if(thisimage.ne.1)then
!   j=startrow
!   k=endrow
!   x(j:k)[1]=x(j:k)[1]+x(j:k)
! endif

 sync all     !needed to wait for the update   ==>  probably better to use sync images

 write(this%unlog,'(/a,i0,a)')' **End image ',thisimage,' of the PCG Coarray subroutine...'
 !$ write(this%unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2

end subroutine

!*FINAL
subroutine destroy_fpcg(this)
 type(fpcg),intent(inout)::this

 call destroy_solver(this)

 this%smalleigenval=0._real64
 this%largeeigenval=0._real64

end subroutine


!PRIVATE
!COEFFICENT MATRIX
subroutine multbyv(this,p,w,startrow,endrow,unlog)
 class(*),intent(inout)::this
 integer(kind=int32),intent(in)::unlog
 integer(kind=int32),intent(in)::startrow,endrow
 real(kind=real64),intent(in)::p(:)
 real(kind=real64),intent(inout)::w(:)
 

 select type(this)
  class is(gen_coeff)
   call this%multbyv(p,w,startrow,endrow)  
  type is(crssparse)
   call this%multbyv(1._real64,'n',p,0._real64,w)
  class default
   write(unlog,'(a)')' ERROR: the proposed type(class) of preconditioner is not supported!'
   error stop
 end select

end subroutine

!CONDITION NUMBER
subroutine addalphabetatot(T,iter,previousalpha,beta)
 integer(kind=int32),intent(in)::iter
 real(kind=real64),intent(in)::previousalpha,beta
 real(kind=real64),intent(inout)::T(:,:)
 
 integer(kind=int32)::previousiter
 real(kind=real64)::b

 b=sqrt(beta)
 previousiter=iter-1

 T(previousiter,previousiter)=T(previousiter,previousiter)+1.d0/previousalpha
 T(previousiter,iter)=T(previousiter,iter)+b/previousalpha
 T(iter,previousiter)=T(iter,previousiter)+b/previousalpha
 T(iter,iter)=T(iter,iter)+beta/previousalpha

end subroutine

subroutine eigenvalandcondnumber(T,iter,small,large,unlog)
 integer(kind=int32),intent(in)::iter
 integer(kind=int32),intent(in),optional::unlog
 real(kind=real64),intent(out),optional::small,large
 real(kind=real64),intent(inout)::T(:,:)

 integer(kind=int32)::un
 real(kind=real64),allocatable::eigenvalue(:),ttmp(:,:)

 un=6
 if(present(unlog))un=unlog

 allocate(eigenvalue(iter-2))

 !call eigensyevd(T(2:iter-1,2:iter-1),iter-2,eigenvalue)

 allocate(ttmp,source=T(2:iter-1,2:iter-1))
write(un,*)iter,'aaa',ttmp(iter-2,iter-2)
 call eigensyevd(ttmp,iter-2,eigenvalue)
 deallocate(ttmp)

 if(present(small))small=minval(eigenvalue)
 if(present(large))large=minval(eigenvalue)

 write(un,'(/"  Mininum | maximum eigenvalues        : ",g0.5,a,g0.5)')minval(eigenvalue),' | ',maxval(eigenvalue)
 write(un,'( "   Effective spectral condition number : ",g0.5)')maxval(eigenvalue)/minval(eigenvalue)
 deallocate(eigenvalue)

end subroutine

subroutine eigensyevd(mat,n,w,leigvectors)
 integer(kind=int32),intent(in)::n
 real(kind=real64),intent(inout)::mat(:,:),w(:)
 logical,intent(in),optional::leigvectors

 integer(kind=int32)::info,i,j,lwork,liwork
 integer(kind=int32),allocatable::iwork(:)
 real(kind=real64),allocatable::work(:)
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

!DESTROY
subroutine destroy_solver(this)
 class(solver),intent(inout)::this

 this%neq=-9
 this%maxit=10000
 this%tol=1.e-6
 this%unlog=6

end subroutine

!PRECONDITIONER
subroutine solveprecond(precond,z,r,unlog)
 !solve precond*z=r
 class(*),intent(inout)::precond
 integer(kind=int32),intent(in)::unlog
 real(kind=real64),intent(out)::z(:)
 real(kind=real64),intent(inout)::r(:)

 select type(precond)
  class is(gen_precond)
   call precond%solve(z,r)
  type is(crssparse)
   call precond%solve(z,r)
  class default
   write(unlog,'(a)')' ERROR: the proposed type(class) of preconditioner is not supported!'
   error stop
 end select

end subroutine

!NORM
function norm(vector,starteq,endeq)
 real(kind=real64),intent(in)::vector(:)
 integer(kind=int32),intent(in)::starteq,endeq
 real(kind=real64)::norm
 
 integer(kind=int32)::i

 norm=0.d0
 do i=starteq,endeq
  norm=norm+vector(i)**2
 enddo
 norm=sqrt(norm)

end function

!READ RHS
subroutine readrhs(rhs,crhs,startrow,endrow,unlog)
 !stupid to read a vector like that, but it works
 integer(kind=int32),intent(in)::startrow,endrow,unlog
 real(kind=real64),intent(inout)::rhs(:)
 character(len=*),intent(in)::crhs

 integer(kind=int32)::i,j,io,un
 real(kind=real64)::val
 !$ real(kind=real64)::t1

 write(unlog,'(/a)')' Start to read the rhs'
 !$ t1=omp_get_wtime()

 rhs=0.d0
 open(newunit=un,file=trim(crhs),access='stream',action='read',status='old')!,buffered='yes')
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

 !$ write(unlog,'(a,f0.3,a)')' End of reading the rhs (Elapsed time (s): ',omp_get_wtime()-t1,')'

end subroutine

end module
