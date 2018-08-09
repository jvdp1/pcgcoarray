!ifort prepafile.f90 -o prepafile

program  prepafile
 use modkind
 use modsparse
 implicit none
 integer(kind=intc)::io,un,unin,i,j,k,l,m,n,image,numimages=4
 integer(kind=int4)::nia,startrow,startcol,endrow,endcol
 integer(kind=intnel)::nel
 character(len=3)::typesparse='col'
 character(len=10)::cdummy
 real(kind=real8)::val
 logical::lex
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

 !call sparse%printfile(500)
 write(*,'(a/)')' The matrix is read...'


 inquire(file='param.prepafile',exist=lex)
 if(lex)then
  open(newunit=unin,file='param.prepafile',action='read')
  numimages=0
  do 
   read(unin,*,iostat=io)i,j,k
   if(io.ne.0)exit
   if(k.ne.sparse%get_dimension_1())typesparse='row'
   numimages=numimages+1
  enddo
  rewind(unin)
 endif

 open(newunit=un,file='param.pcgcoarray.'//adjustl(typesparse(:len_trim(typesparse))),status='replace',action='write')
 do i=1,numimages
  if(lex)then
   read(unin,*)image,startrow,endrow,startcol,endcol
  else
   image=i
   startrow=1
   endrow=sparse%get_dimension_1()
   nel=int(real(endrow)/numimages)
   startcol=nel*(image-1)+1
   endcol=startcol+nel-1
   if(image.eq.numimages)endcol=sparse%get_dimension_2()
  endif

  sparsesub=sparse%sub(startrow,endrow,startcol,endcol)

  !call sparsesub%printfile(500+image)
  call sparsesub%printbin(image,typesparse)
  call sparsesub%reset
  write(un,*)image,startrow,endrow,startcol,endcol
 enddo

 close(un)
 if(lex)close(unin)

end program
