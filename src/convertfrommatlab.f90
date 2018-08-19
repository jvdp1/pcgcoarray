program convertfrommatlab
 use modkind
 use modsparse
 implicit none
 
 integer(kind=int4)::i,j,k,unin,unout,un,io,nmat,cim,ndiag
 integer(kind=int4)::dim1,dim2
 integer(kind=int4),allocatable::ia(:)
 real(kind=real8),allocatable::ra(:)
 character(len=30)::cdim,cia,cja,csa
 character(len=30)::cdummy1,cdummy2
 type(csr)::sparse,precond

 !per image
 open(newunit=unin,file='param.convert',action='read',status='old')

 call numlines(unin,nmat)
 nmat=nmat-1

 open(newunit=unout,file='param.pcgcoarray.row',status='replace',action='write')

 do i=1,nmat
  read(unin,*,iostat=io)cim,cdim,cia,cja,csa
  if(io.ne.0)exit
  write(*,'(/a,i0)')' Image ',cim
  !Dimension
  open(newunit=un,file=cdim,action='read',status='old')
  do j=1,3
   read(un,*)
  enddo
  read(un,*)dim1
  do j=1,4
   read(un,*)
  enddo
  read(un,*)dim2
  close(un)
  write(*,*)'dim1 x dim2: ',dim1,dim2
  !read ia
  open(newunit=un,file=cia,action='read',status='old')
  do j=1,3
   read(un,*)
  enddo
  read(un,*)cdummy1,cdummy2,k
  allocate(ia(k))
  read(un,*)
  do j=1,k
   read(un,*)ia(j)
  enddo
  close(un)
  !allocate sparse
  call sparse%alloc(ia(k)-1,dim1,dim2)
  sparse%ia=ia
  deallocate(ia)
  !read ja
  open(newunit=un,file=cja,action='read',status='old')
  do j=1,3
   read(un,*)
  enddo
  read(un,*)cdummy1,cdummy2,k
  allocate(ia(k))
  read(un,*)
  do j=1,k
   read(un,*)ia(j)
  enddo
  close(un)
  sparse%ja=ia
  deallocate(ia)
  !read ra
  open(newunit=un,file=csa,action='read',status='old')
  do j=1,3
   read(un,*)
  enddo
  read(un,*)cdummy1,cdummy2,k
  allocate(ra(k))
  read(un,*)
  do j=1,k
   read(un,*)ra(j)
  enddo
  close(un)
  sparse%a=ra
  deallocate(ra)
  ndiag=0 !dim2
  precond=sparse%subtriu(1,dim1,(i-1)*dim1+1,i*dim1,ndiag+(i-1)*dim1,(i-1)*dim1)
  call precond%printbin(cim,'precond','row')
  call sparse%printbin(cim,'row')
  call sparse%reset()
  call precond%reset()
  write(unout,*)cim,(i-1)*dim1+1,i*dim1,1,dim2
 enddo
 close(unout)
 close(unin)

 allocate(ra(dim2))
 open(newunit=un,file='rhs.dat',action='read',status='old')
 do j=1,5
  read(un,*)
 enddo
 do j=1,dim2
  read(un,*)ra(j)
 enddo
 close(un)
 open(newunit=un,file='rhs.bin',access='stream',action='write',status='replace',buffered='yes')
 write(un)dim2
 write(un)ra
 close(un)

contains

subroutine numlines(unfile,n)
 integer,intent(in)::unfile
 integer,intent(out)::n

 integer::io
 rewind(unfile)
 n=0
 do
  read(unfile,*,iostat=io)
  if (io.ne.0) exit
  n=n+1
 enddo
 rewind(unfile)
end subroutine

end program
