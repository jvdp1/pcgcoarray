!ifort prepafile.f90 -o prepafile

program  prepafile
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::io,un,i,j,k,l,m,n,numimages=4
 integer(kind=int4)::nia,startrow,startcol,endrow,endcol
 integer(kind=intnel)::nel
 character(len=10)::cdummy
 real(kind=real8)::val
 type(csr)::sparse,sparsesub

 open(newunit=un,file='matrixija.bin',status='old',action='read',access='stream')
 read(un)nia
 read(un)nel
 call sparse%alloc(nel,nia-1)
 
 read(un)sparse%ia
 print*,'read ia'
 read(un)sparse%ja
 print*,'read ja'
 read(un)sparse%a
 print*,'read a'
 close(un)

 call sparse%printfile(500)
 write(*,'(a/)')' The matrix is read...'

 open(newunit=un,file='param.pcgcoarray',status='replace',action='write')
 do i=1,numimages
  startrow=1
  endrow=sparse%n
  nel=int(real(endrow)/numimages)
  startcol=nel*(i-1)+1
  endcol=startcol+nel-1
  if(i.eq.numimages)endcol=sparse%m

  sparsesub=sparse%sub(startrow,endrow,startcol,endcol)

!  call sparsesub%printfile(500+i)
  call sparsesub%printbin(i)
  call sparsesub%reset
  write(un,*)i,startrow,endrow,startcol,endcol
 enddo
 close(un)



end program
