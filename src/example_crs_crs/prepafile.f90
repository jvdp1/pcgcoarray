program  prepafile
 use iso_fortran_env
 use modsparse
 implicit none
 integer(kind=int32)::io,un,unin,i,j,k,l,m,n,image,numimages=4
 integer(kind=int32)::nia,startrow,startcol,endrow,endcol,ndiag
 integer(kind=int32)::nel
 integer(kind=int32),allocatable::ia(:),ja(:)
 character(len=3)::typesparse='col'
 character(len=100)::cdummy
 real(kind=real64)::val
 real(kind=real64),allocatable::a(:)
 logical::lex
 type(crssparse)::crs,crssub,crssubtmp

 open(newunit=un,file='matrixija.bin',status='old',action='read',access='stream')
 read(un)nia
 read(un)nel
 
 allocate(ia(nia),ja(nel),a(nel))

 crs=crssparse(nia-1,nel,lupper=.false.)

 read(un)ia
 print*,'read ia'
 read(un)ja
 print*,'read ja'
 read(un)a
 print*,'read a'
 close(un)

 call crs%external(ia,ja,a)
 deallocate(ia,ja,a)

 call crs%printstats()

 write(*,'(a/)')' The matrix is read...'


 inquire(file='param.prepafile',exist=lex)
 if(lex)then
  open(newunit=unin,file='param.prepafile',action='read')
  numimages=0
  do 
   read(unin,*,iostat=io)i,j,k
   if(io.ne.0)exit
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
   endrow=crs%getdim(1)
   nel=int(real(endrow)/numimages)
   startcol=nel*(image-1)+1
   endcol=startcol+nel-1
   if(image.eq.numimages)endcol=crs%getdim(2)
   ndiag=0
  endif

  crssub=crs%submatrix(startrow,endrow,startcol,endcol)
  call crssub%printstats()

  write(un,*)image,startrow,endrow,startcol,endcol
 
  
  write(cdummy,'(i0)')image
  cdummy='crs_subpcg.'//adjustl(typesparse(:len_trim(typesparse)))//adjustl(cdummy(:len_trim(cdummy)))
  call crssub%save(cdummy)
  call crssub%destroy()


  !preconditioner
  if((endrow-startrow).lt.endcol-startcol)then     !row format

   crssubtmp=crs%submatrix(startrow,endrow,startrow,endrow)
   crssub=crssubtmp%diag(ndiag)
   call crssubtmp%destroy()
   call crssub%sort()
!   write(cdummy,'(i0)')image
!   cdummy='txt_precond.'//adjustl(typesparse(:len_trim(typesparse)))//adjustl(cdummy(:len_trim(cdummy)))
!   call crssub%printtofile(cdummy)
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
