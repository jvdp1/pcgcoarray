program  pcgcorray
 !$ use omp_lib
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::thisimage
 integer(kind=intc)::i,j,k,l,m,n,neq,ncol,iter
 integer(kind=int4)::startcolk,maxit=4000
 integer(kind=int4)::startrow[*],endrow[*],startcol[*],endcol[*]
 integer(kind=intnel)::nel
 character(len=80)::host,cdummy
 real(kind=real8)::tol=1.e-12
 real(kind=real8)::oldtau,conv,thr,beta
 real(kind=real8)::b_norm[*],resvec1[*],alpha[*],tau[*]
 real(kind=real8),allocatable::rhs(:),precond(:)
 real(kind=real8),allocatable::z(:),p(:)
 real(kind=real8),allocatable::x(:)[:],r(:)[:],w(:)[:]
 real(kind=real8)::t1,val
 type(csr)::sparse

 thisimage=this_image()

 write(cdummy,'(i0)')thisimage
 if(thisimage.eq.1)then
  write(*,'(/a/)')' PCG solver with a LHS divided by columns...'
 endif

 call get_environment_variable("HOSTNAME",value=host)
 write(*,'(2(a,i0),2a)')"Hello from image ",thisimage," out of ",num_images()," total images on host ",trim(host)

 !read the parameter file on image 1
 if(thisimage.eq.1)then
  call readparam(startrow,endrow,startcol,endcol)
 endif

 sync all

 !Reads the matrix
 call sparse%get('subpcg.col'//adjustl(cdummy(:len_trim(cdummy))))

 !call sparse%printfile(600+thisimage)

 !sync all

 !read rhs
 allocate(rhs(sparse%get_dimension_2()))
 call readrhs(rhs,startcol,endcol)

 write(*,'(a,i0)')' Preparation for the PCG for image ',thisimage
 neq=sparse%get_dimension_1()
 ncol=sparse%get_dimension_2()

 !create preconditioner
 precond=sparse%diagcol(startcol)

 !invert preconditioner
 do i=1,ncol
  if(precond(i).ne.0)precond(i)=1.d0/precond(i)
 enddo
 
 !Initialistion PCG arrays
 allocate(x(neq)[*])
 x=0.d0

 allocate(z(ncol),p(ncol))
 allocate(r(neq)[*],w(neq)[*])
 r=0.d0;p=0.d0;z=0.d0;w=0.d0

 oldtau=1.d0

 !initiatilisation of r
 if(sum(x).eq.0.d0)then
  r(startcol:endcol)=rhs
 else
  print*,'not yet implemented'
  error stop
 endif

 b_norm=norm(rhs,1,ncol)
 resvec1=norm(r,startcol,endcol)

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
 !sync all  !not sure if it is really needed
 do i=1,num_images()
  if(i.ne.thisimage)then
   !receives updates from other for its own image
   j=startcol[i]
   k=endcol[i]
   r(j:k)=r(j:k)+r(j:k)[i]
  endif
 enddo
 sync all  !not sure if it is really needed

 conv=resvec1/b_norm

 if(thisimage.eq.1)then
  write(*,'(a,e15.5)')' Norm of RHS: ',b_norm
  write(*,'(" Iteration ",i6," Convergence = ",e12.5)')1,conv
 endif

 alpha=1.d0
 iter=2
 thr=tol*b_norm

 !Start iteration
 !$ t1=omp_get_wtime() 
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
  do i=1,ncol
   p(i)=z(i)+beta*p(i)
  enddo

  !w=LHS*p
  call sparse%multv(p,w)
 
  !update w
  sync all
  if(thisimage.eq.1)then
   do i=2,num_images()
    !receives updates from other for its own image
   w=w+w(:)[i]
   enddo
  endif
  sync all
  !2. update on the other image
  if(thisimage.ne.1)w(:)=w(:)[1]

  !alpha=p*w
  alpha=0.d0
  do i=1,ncol
   alpha=alpha+p(i)*w(startcolk+i)
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

  do i=1,ncol
   x(startcolk+i)=x(startcolk+i)+alpha*p(i)
  enddo

  do i=1,neq
   r(i)=r(i)-alpha*w(i)
  enddo
  resvec1=norm(r,1,neq)

  conv=resvec1/b_norm


  if(thisimage.eq.1)then
   write(*,'(" Iteration ",i6," Convergence = ",e12.5,x,e12.5)')iter,conv,resvec1
  endif

  iter=iter+1

 enddo
 !$ if(thisimage.eq.1)then
 !$ val=omp_get_wtime()-t1
 !$  write(*,'("  Wall clock time for the iterative process (seconds): ",f12.2)')val
 !$  write(*,'("  Approximate Wall clock time per iteration (seconds): ",f12.2)')val/(iter-1)
 !$ endif

 sync all

 if(thisimage.eq.1)then
  do i=2,num_images()
   !receives updates from other for its own image
   j=startcol[i]
   k=endcol[i]
   x(j:k)=x(j:k)+x(j:k)[i]
  enddo
 endif
 sync all
 if(thisimage.eq.1)call print_ascii(x,1,neq)


 write(*,'(2(a,i0),a)')"End for image ",thisimage," out of ",num_images()," total images!"

contains

subroutine readparam(startrow,endrow,startcol,endcol)
 integer(kind=int4),intent(inout)::startrow[*],endrow[*],startcol[*],endcol[*]

 integer(kind=int4)::io,un,i,j,k,l,m,n

 open(newunit=un,file='param.pcgcoarray.col',status='old',action='read')
 n=0
 do
  read(un,*,iostat=io)i,j,k,l,m
  if(io.ne.0)exit
  n=n+1
  if(n.gt.num_images())then
   write(*,'(a)')' The parameter file is not correct!'
   error stop
  endif
  startrow[i]=j
  endrow[i]=k
  startcol[i]=l
  endcol[i]=m
 enddo
 close(un)
 if(n.ne.num_images())then
  write(*,'(a)')' The parameter file is not correct!'
  error stop
 endif

end subroutine

subroutine readrhs(rhs,startcol,endcol)
 !stupid to read a vector like that, but it works
 integer(kind=int4),intent(in)::startcol,endcol
 real(kind=real8),intent(inout)::rhs(:)

 integer(kind=int4)::i,j,io,un
 real(kind=real8)::val

 write(*,'(/a)')' Starts to read the rhs'
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
 write(*,'(a)')' End of reading the rhs'

end subroutine

subroutine print_ascii(x,startpos,endpos)
 integer(kind=intc),intent(in)::startpos,endpos
 real(kind=real8),intent(in)::x(:)

 integer(kind=intc)::un

 write(*,'(a,i0,a)')" Image ",this_image()," starts to write the solutions!"

#if (TOUTPUT==1)
 open(newunit=un,file='solutions.pcgcoarray.col.dist')
#elif (TOUTPUT==2)
 open(newunit=un,file='solutions.pcgcoarray.col.shared')
#else
 open(newunit=un,file='solutions.pcgcoarray.col')
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
