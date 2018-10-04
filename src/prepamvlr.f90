program prepamvlr
 use modkind
 use modmvlr
 implicit none
 integer(kind=int4)::i,j,io,un
 integer(kind=int4)::ntrait,neffect,ncov,nphen
 integer(kind=int4)::unlog=6
 integer(kind=int4)::unout
 integer(kind=int4)::neq
 integer(kind=int4),allocatable::itmp(:),maxlevel(:),startlevel(:),model(:,:)
 integer(kind=real8)::itmp8
 real(kind=real8)::missing
 real(kind=real8),allocatable::val(:)
 real(kind=real8),allocatable::array(:)
 character(len=20)::datafile,dataoutput='data.eq'
 type(mvlr)::reg

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

 !sum of levels
 allocate(startlevel(neffect))
 startlevel=0
 do i=2,neffect
  startlevel(i)=startlevel(i-1)+maxlevel(i-1)
 enddo
 startlevel=startlevel*ntrait
 print*,'startlevel ',startlevel

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

 !compute preconditioner
 call reg%init('parammvlr.dat',500+this_image())
 array=reg%diag()
 neq=reg%getneq()

 open(newunit=unout,file='precond.stream',access='stream',status='replace',action='write')
 write(unout)neq
 write(unout)array
 close(unout)

contains

function getaddress(startlevel,ntrait,effect,level,trait) result(address)
 integer(kind=4)::address
 integer(kind=4),intent(in)::startlevel(:),ntrait,effect,level,trait
 
 address=startlevel(effect)+(level-1)*ntrait+trait

end function

function setmissing(val,missing,ntrait) result(combi)
 integer(kind=int4),intent(in)::ntrait
 integer(kind=real8)::combi
 real(kind=real8),intent(in)::val(:),missing

 integer(kind=int4)::i
 
 combi=0
 do i=1,ntrait
  if(val(i).ne.missing)combi=ibset(combi,i-1)
 enddo

end function

end program
