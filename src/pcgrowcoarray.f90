program  pcgrowcorray
 !$ use omp_lib
 !$ use mkl_service
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::thisimage,unlog,unconv
 integer(kind=intc)::i,j,k,l,m,n,neq,nrow,iter
 integer(kind=int4)::startrowk,maxit=4000
 integer(kind=int4)::startrow[*],endrow[*],startcol[*],endcol[*]
 integer(kind=intnel)::nel
 character(len=80)::host,cdummy
 real(kind=real8)::tol=1.e-12
 real(kind=real8)::oldtau,conv,thr,beta
 real(kind=real8)::b_norm[*],resvec1[*],alpha[*],tau[*]
 real(kind=real8),allocatable::rhs(:),precond(:)
 real(kind=real8),allocatable::z(:),r(:),w(:)
 real(kind=real8),allocatable::x(:)[:],p(:)[:]
 !$ real(kind=real8)::t1,t2,val
 type(csr)::sparse

 !$ t2=omp_get_wtime() 

 thisimage=this_image()

 write(cdummy,'(i0)')thisimage
 open(newunit=unlog,file='log.pcgrow.'//adjustl(cdummy(:len_trim(cdummy))),action='write',status='replace')

 write(unlog,'(/a)')' PCG solver with a LHS divided by rows...'

 call get_environment_variable("HOSTNAME",value=host)
 write(unlog,'(/2(a,i0),2a)')" Hello from image ",thisimage," out of ",num_images()," total images on host ",trim(host)

 write(unlog,'(/" Number of images            : ",i0)')num_images()
 !$omp parallel
 !$omp master
 !$ write(unlog,'(" Number of threads for OpenMP: ",i0)')omp_get_num_threads() 
 !$omp end master
 !$omp end parallel
 !$ write(unlog,'(" Number of threads for MKL   : ",i0)')mkl_get_max_threads() 

 !read the parameter file on image 1
 if(thisimage.eq.1)then
  call readparam(startrow,endrow,startcol,endcol,unlog)
 endif

 sync all

 !Reads the matrix
 call sparse%get('subpcg.row'//adjustl(cdummy(:len_trim(cdummy))),unlog)

 !call sparse%printfile(600+thisimage)

 !sync all

 !read rhs
 allocate(rhs(sparse%get_dimension_1()))
 call readrhs(rhs,startrow,endrow,unlog)

 write(unlog,'(/a,i0)')' Preparation for the PCG for image ',thisimage
 neq=sparse%get_dimension_2()
 nrow=sparse%get_dimension_1()

 !create preconditioner
 precond=sparse%diagrow(startrow)

 !invert preconditioner
 do i=1,nrow
  if(precond(i).ne.0)precond(i)=1.d0/precond(i)
 enddo
 
 !Initialistion PCG arrays
 allocate(x(neq)[*])
 x=0.d0

 allocate(z(nrow),w(nrow),r(nrow))
 allocate(p(neq)[*])
 r=0.d0;p=0.d0;z=0.d0;w=0.d0

 oldtau=1.d0

 !initiatilisation of r
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
  write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,x,e12.5)')1,conv,0.
 endif

 alpha=1.d0
 iter=2
 thr=tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
 startrowk=startrow-1
 do while(resvec1.gt.thr.and.iter.le.maxit)
  !z=M*r
  do i=1,nrow
   z(i)=precond(i)*r(i)
  enddo

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
  call sparse%multv(p,w)
 
  !alpha=p*w
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
   write(unconv,'(" Iteration ",i6," Convergence = ",e12.5,x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 !$ if(thisimage.eq.1)then
 !$  val=omp_get_wtime()-t1
 !$  write(unlog,'(/"  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(unlog,'("  Approximate wall clock time per iteration (seconds): ",f12.2)')val/(iter-1)
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
 if(thisimage.eq.1)call print_ascii(x,1,neq,unlog)

 write(unlog,'(/2(a,i0),a)')" End for image ",thisimage," out of ",num_images()," total images!"
 !$ write(unlog,'("   Wall clock time: ",f12.2)')omp_get_wtime()-t2
 close(unlog)

contains

subroutine readparam(startrow,endrow,startcol,endcol,unlog)
 integer(kind=int4),intent(in)::unlog
 integer(kind=int4),intent(inout)::startrow[*],endrow[*],startcol[*],endcol[*]

 integer(kind=int4)::io,un,i,j,k,l,m,n

 open(newunit=un,file='param.pcgcoarray.row',status='old',action='read')
 n=0
 do
  read(un,*,iostat=io)i,j,k,l,m
  if(io.ne.0)exit
  n=n+1
  if(n.gt.num_images())then
   write(unlog,'(a)')' The parameter file is not correct!'
   error stop
  endif
  startrow[i]=j
  endrow[i]=k
  startcol[i]=l
  endcol[i]=m
 enddo
 close(un)
 if(n.ne.num_images())then
  write(unlog,'(a)')' The parameter file is not correct!'
  error stop
 endif

end subroutine

subroutine readrhs(rhs,startrow,endrow,unlog)
 !stupid to read a vector like that, but it works
 integer(kind=int4),intent(in)::startrow,endrow,unlog
 real(kind=real8),intent(inout)::rhs(:)

 integer(kind=int4)::i,j,io,un
 real(kind=real8)::val

 write(unlog,'(/a)')' Starts to read the rhs'
 rhs=0.d0
 open(newunit=un,file='rhs.bin',access='stream',action='read',status='old',buffered='yes')
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

subroutine print_ascii(x,startpos,endpos,unlog)
 integer(kind=intc),intent(in)::startpos,endpos,unlog
 real(kind=real8),intent(in)::x(:)

 integer(kind=intc)::un

 write(unlog,'(/a,i0,a)')" Image ",this_image()," starts to write the solutions!"

#if (TOUTPUT==1)
 open(newunit=un,file='solutions.pcgcoarray.row.dist')
#elif (TOUTPUT==2)
 open(newunit=un,file='solutions.pcgcoarray.row.shared')
#else
 open(newunit=un,file='solutions.pcgcoarray.row')
#endif

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
