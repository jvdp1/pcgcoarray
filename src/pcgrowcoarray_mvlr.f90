program  pcgrowcorray_mvlr
 !$ use omp_lib
 !$ use mkl_service
 use modkind
 use modsparse
 use modpcgcoarray
 use modprecond
 use modmvlr
 implicit none
 integer(kind=int4)::thisimage,unlog
 integer(kind=int4)::neq
 integer(kind=int4)::un
 integer(kind=int4)::startrow[*],endrow[*],startcol[*],endcol[*]
 character(len=80)::host,cdummy
 real(kind=real8),allocatable::array(:)
 real(kind=real8),allocatable::x(:)[:]
 !$ real(kind=real8)::t2
 type(arrayprecond)::precond
 type(mvlr)::reg

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
 cdummy='parammvlr.dat'
 call reg%init(cdummy,unlog)

 neq=reg%getneq()


 !create preconditioner
 write(unlog,'(/a)')' Extraction of the diagonal elements...'

 precond%dim1=endrow-startrow+1
 allocate(precond%array(precond%dim1))
 open(newunit=un,file='precond.stream',access='stream',status='old',action='read')
 read(un)neq
 allocate(array(neq))
 read(un)array
 close(un)
 precond%array=array(startrow:endrow)
 deallocate(array)

 !solution vector
 allocate(x(neq)[*])
 x=0.d0

 sync all

 call pcgrowcoarray(neq,reg,x,'rhs.bin',precond,startrow,endrow,unlog)

 sync all

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

!PRINT
subroutine print_ascii(x,startpos,endpos,unlog)
 integer(kind=int4),intent(in)::startpos,endpos,unlog
 real(kind=real8),intent(in)::x(:)

 integer(kind=int4)::un,i

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

end program
