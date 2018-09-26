program  prepafile
 use modkind
 use modsparse
 !use modsparse_old
 implicit none
 integer(kind=intc)::io,un,unin,i,j,k,l,m,n,image,numimages=4
 integer(kind=int4)::nia,startrow,startcol,endrow,endcol,ndiag
 integer(kind=int4)::nel
 integer(kind=int4),allocatable::ia(:),ja(:)
 character(len=3)::typesparse='col'
 character(len=100)::cdummy
 real(kind=real8)::val
 real(kind=real8),allocatable::a(:)
 logical::lex
! type(csr)::sparse,sparsesub
 type(crssparse)::crs,crssub,crssubtmp

 open(newunit=un,file='matrixija.bin',status='old',action='read',access='stream')
 read(un)nia
 read(un)nel
! call sparse%alloc(nel,nia-1)
 
 allocate(ia(nia),ja(nel),a(nel))

 crs=crssparse(nia-1,nel,lupper=.false.)

 read(un)ia
! sparse%ia=ia
 print*,'read ia'
 read(un)ja
! sparse%ja=ja
 print*,'read ja'
 read(un)a
! sparse%a=a
 print*,'read a'
 close(un)

 call crs%external(ia,ja,a)
 deallocate(ia,ja,a)

 call crs%printstats()

 !call sparse%printfile(500)
 write(*,'(a/)')' The matrix is read...'


 inquire(file='param.prepafile',exist=lex)
 if(lex)then
  open(newunit=unin,file='param.prepafile',action='read')
  numimages=0
  do 
   read(unin,*,iostat=io)i,j,k
   if(io.ne.0)exit
   !if(k.ne.sparse%get_dimension_1())typesparse='row'
   if(k.ne.crs%getdim(1))typesparse='row'
   numimages=numimages+1
  enddo
  rewind(unin)
 endif

 open(newunit=un,file='param.pcgcoarray.'//adjustl(typesparse(:len_trim(typesparse))),status='replace',action='write')
 do i=1,numimages
  if(lex)then
   read(unin,*)image,startrow,endrow,startcol,endcol,ndiag
  else
   image=i
   startrow=1
   !endrow=sparse%get_dimension_1()
   endrow=crs%getdim(1)
   nel=int(real(endrow)/numimages)
   startcol=nel*(image-1)+1
   endcol=startcol+nel-1
   !if(image.eq.numimages)endcol=sparse%get_dimension_2()
   if(image.eq.numimages)endcol=crs%getdim(2)
   ndiag=0
  endif

!  sparsesub=sparse%sub(startrow,endrow,startcol,endcol)
  crssub=crs%submatrix(startrow,endrow,startcol,endcol)
  call crssub%printstats()

  !call sparsesub%printfile(500+image)
!  call sparsesub%printbin(image,typesparse)
!  call sparsesub%reset
  write(un,*)image,startrow,endrow,startcol,endcol
 
  
  write(cdummy,'(i0)')image
  cdummy='crs_subpcg.'//adjustl(typesparse(:len_trim(typesparse)))//adjustl(cdummy(:len_trim(cdummy)))
  call crssub%save(cdummy)
  call crssub%destroy()


  !preconditioner
  if((endrow-startrow).lt.endcol-startcol)then     !row format
!   !sparsesub=sparse%subdiag(startrow,endrow,startrow,endrow)
!   !sparsesub=sparse%subup(startrow,endrow,startrow,endrow)
!   sparsesub=sparse%subtriu(startrow,endrow,startrow,endrow,ndiag)
!   call sparsesub%sort()
!   !call sparsesub%printfile(650+image)
!   call sparsesub%printbin(image,'precond',typesparse)
!   call sparsesub%reset

   crssubtmp=crs%submatrix(startrow,endrow,startrow,endrow)
   crssub=crssubtmp%diag(ndiag)
   call crssubtmp%destroy()
   call crssub%sort()
   write(cdummy,'(i0)')image
   cdummy='crs_precond.'//adjustl(typesparse(:len_trim(typesparse)))//adjustl(cdummy(:len_trim(cdummy)))
   call crssub%printstats()
   call crssub%save(cdummy)
   call crssub%destroy()
  elseif((endrow-startrow).gt.endcol-startcol)then !column format
   write(*,'(a)')' ERROR: preconditioner option not supported'
   stop
   !sparsesub=sparse%subup(startcol,endcol,startcol,endcol)
  endif

 enddo

 close(un)
 if(lex)close(unin)

end program
