program prepamvlr
 use iso_fortran_env
 use modmvlr
 use modsparse
 implicit none
 integer(kind=int32)::i,j,io,un
 integer(kind=int32)::ntrait,neffect,ncov,nphen
 integer(kind=int32)::unlog=6
 integer(kind=int32)::unout
 integer(kind=int32)::neq
 integer(kind=int32)::nsparse
 integer(kind=int32),allocatable::itmp(:),maxlevel(:),startlevel(:),model(:,:)
 integer(kind=real64)::itmp8
 real(kind=real64)::missing
 real(kind=real64),allocatable::val(:)
 real(kind=real64),allocatable::array(:)
 character(len=20)::datafile,dataoutput='data.eq'
 character(len=20),allocatable::namesparse(:)
 type(mvlr)::reg
 type(crssparse)::mat

 write(unlog,'(a)')' Name of data file?'
 read(*,*)datafile
 write(unlog,'(a)')' Number of traits?'
 read(*,*)ntrait
 write(unlog,'(a)')' Number of effects?'
 read(*,*)neffect
 write(unlog,'(a)')' Number of covariables?'
 read(*,*)ncov
 write(unlog,'(a)')' Missing phenotype?'
 read(*,*)missing

 !Number of levels per effect
 allocate(itmp(neffect),maxlevel(neffect),val(ntrait))
 maxlevel=-1
 nphen=0
 open(newunit=un,file=datafile,action='read')
 do
  read(un,*,iostat=io)val,itmp
  if(io.ne.0)exit
  do i=1,neffect
   if(itmp(i).gt.maxlevel(i))maxlevel(i)=itmp(i)
  enddo
  nphen=nphen+1
 enddo

 !Model per trait
 allocate(model(neffect,ntrait))
 do i=1,ntrait
  write(unlog,'(a,i0,a)')' Model for trait ',i,'?'
  read(*,*,iostat=io)model(:,i)
  if(io.ne.0)then
   write(unlog,*)' ERROR'
   stop
  endif
 enddo

 write(unlog,*)' Max levels: ',maxlevel
 write(unlog,*)' # pheno   : ',nphen

 !possible additions of CRS matrices
 write(unlog,'(a)')' Number of sparse matrices?'
 read(*,*)nsparse
 
 if(nsparse.gt.0)then 
  allocate(namesparse(nsparse))
  namesparse=''
  do i=1,nsparse
   write(unlog,'(a,i0,a)')' Name of sparse matrix file ',i,'?'
   read(*,*)namesparse(i)
  enddo
 endif


 !sum of levels
 allocate(startlevel(neffect))
 startlevel=0
 do i=2,neffect
  startlevel(i)=startlevel(i-1)+maxlevel(i-1)
 enddo
 startlevel=startlevel*ntrait
 write(unlog,*)'startlevel ',startlevel

 rewind(un)
 open(newunit=unout,file=dataoutput,action='write',status='replace')
 do
  read(un,*,iostat=io)val,itmp
  if(io.ne.0)exit
  do i=1,neffect
   itmp(i)=getaddress(startlevel,ntrait,i,itmp(i),1)
  enddo
  itmp8=setmissing(val,missing,ntrait)
  write(unout,'(<ntrait>(g0.4,1x),i0,1x,<neffect>(i0,1x))')val,itmp8,itmp
 enddo
 close(unout)
 deallocate(val,itmp)

 open(newunit=unout,file='parammvlr.dat',action='write',status='replace')
 write(unout,'(a,t20,a)')trim(dataoutput),'datafile'
 write(unout,'(i0,t20,a)')nphen,'nphen'
 write(unout,'(i0,t20,a)')ntrait,'ntrait'
 write(unout,'(i0,t20,a)')neffect,'neffect'
 write(unout,'(i0,t20,a)')ncov,'ncov'
 write(unout,'(i0,t20,a)')getaddress(startlevel,ntrait,neffect,maxlevel(neffect),ntrait),'neq'
 write(unout,'(a,t20,a)')'res.dat','residual variance file'
 do i=1,ntrait
  write(unout,'(<neffect>i2,t30,a,i0)')model(:,i),'model trait ',i
 enddo
 close(unout)

 open(newunit=unout,file='paramsparse.dat',action='write',status='replace')
 write(unout,'(i0,t20,a)')nsparse,'datafile'
 do i=1,nsparse
  write(unout,'(a,t20,a,i0)')namesparse(i),'sparse ',i
 enddo
 close(unout)


 !compute preconditioner
 call reg%init('parammvlr.dat',500+this_image())
 array=reg%diag()
 neq=reg%getneq()

 !add diagonal elements of sparse matrices
 if(nsparse.gt.0)then
  write(unlog,'(a)')' Addition of diagonal elements of sparse matrices'
  do i=1,nsparse
   mat=crssparse(namesparse(i))
   call mat%printstats()
   array=array+mat%diag()
  enddo
 endif


 open(newunit=unout,file='precond.stream',access='stream',status='replace',action='write')
 write(unout)neq
 write(unout)array
 close(unout)

 close(unlog)

contains

function getaddress(startlevel,ntrait,effect,level,trait) result(address)
 integer(kind=4)::address
 integer(kind=4),intent(in)::startlevel(:),ntrait,effect,level,trait
 
 address=startlevel(effect)+(level-1)*ntrait+trait

end function

function setmissing(val,missing,ntrait) result(combi)
 integer(kind=int32),intent(in)::ntrait
 integer(kind=real64)::combi
 real(kind=real64),intent(in)::val(:),missing

 integer(kind=int32)::i
 
 combi=0
 do i=1,ntrait
  if(val(i).ne.missing)combi=ibset(combi,i-1)
 enddo

end function

end program
