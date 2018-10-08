module modmvlr
 !$ use omp_lib
 use modkind
 use modsparse
#if (COARRAY==1)
 use modcoarraysolver
#endif

 implicit none
 private
 public::mvlr

 type::typeres
  real(kind=real8),allocatable::mat(:,:)
 end type

#if (COARRAY==1)
 type,extends(gen_coeff)::mvlr
#else
 type::mvlr
#endif
  private
  integer(kind=int4)::unlog=6
  integer(kind=int4)::nphen,ntrait,neffect,ncov,neq
  integer(kind=int4)::startphen,endphen,nsubphen
  integer(kind=int4),allocatable::effect(:,:)
  integer(kind=int8),allocatable::presentpheno(:)
  real(kind=real8),allocatable::pheno(:,:)
  real(kind=real8),allocatable::resmati(:,:)
  character(len=30)::datafile,resfile
  logical,allocatable::model(:,:)
  type(typeres),allocatable::tres(:)
 contains
  private
  procedure,public::init=>init_mvlr
  procedure,public::getneq
#if (COARRAY==1)
  procedure,public::multbyv=>multbyv_mvlr
#else
  procedure,public::diag=>computediag_vect_mvlr
#endif
 end type

contains

!PUBLIC
!**GET #EQUATIONS
function getneq(this) result(neq)
 class(mvlr),intent(in)::this
 integer(kind=int4)::neq

 neq=this%neq
end function

!**INIT
subroutine init_mvlr(this,paramfile,unlog)
 class(mvlr),intent(inout)::this
 character(len=*),intent(in)::paramfile
 integer(kind=int4),intent(in),optional::unlog

 integer(kind=int4)::i,k,un,io,nsize
 integer(kind=int4),allocatable::itmp(:)
 integer(kind=int8)::i8tmp,maxcombi
 real(kind=real8),allocatable::rtmp(:)
 logical,allocatable::ltmp(:)

 if(present(unlog))this%unlog=unlog
 
 open(newunit=un,file=paramfile,status='old',action='read')
 read(un,*)this%datafile
 read(un,*)this%nphen
 read(un,*)this%ntrait
 read(un,*)this%neffect
 read(un,*)this%ncov
 read(un,*)this%neq
 read(un,*)this%resfile
 allocate(this%model(this%ntrait,this%neffect),itmp(this%neffect))
 this%model=.false.
 do i=1,this%ntrait
  read(un,*)itmp
  where(itmp.ne.0)this%model(i,:)=.true.
 enddo
 deallocate(itmp)
 close(un)

 write(this%unlog,'(/a,t20,a)')this%datafile,'datafile'
 write(this%unlog,'(i0,t20,a)')this%nphen,'nphen'
 write(this%unlog,'(i0,t20,a)')this%ntrait,'ntrait'
 write(this%unlog,'(i0,t20,a)')this%neffect,'neffect'
 write(this%unlog,'(i0,t20,a)')this%ncov,'ncov'
 write(this%unlog,'(i0,t20,a)')this%neq,'neq'
 write(this%unlog,'(a,t20,a)')this%resfile,'residual variance file'
 allocate(ltmp(this%neffect))
 do i=1,this%ntrait
  ltmp=this%model(i,:)
  write(this%unlog,'(<this%neffect>l1,t30,a,i0)')ltmp,'model trait ',i
 enddo
 deallocate(ltmp)
 
 !start and end data
 nsize=int(real(this%nphen)/num_images())
 this%startphen=nsize*(this_image()-1)+1
 this%endphen=this%startphen+nsize-1
 if(this_image().eq.num_images())this%endphen=this%nphen
 write(this%unlog,'(i0,t20,a)')this%startphen,'start pheno'
 write(this%unlog,'(i0,t20,a)')this%endphen,'end pheno'
 this%nsubphen=this%endphen-this%startphen+1

 !read data file
 
 allocate(this%pheno(this%ntrait,this%nsubphen),this%effect(this%neffect,this%nsubphen)&
           ,this%presentpheno(this%nsubphen))
 allocate(rtmp(this%ntrait),itmp(this%neffect))
 open(newunit=un,file=this%datafile,status='old',action='read')
 k=0
 do i=1,this%nphen
  read(un,*,iostat=io)rtmp,i8tmp,itmp
  if(i.lt.this%startphen)cycle
  if(i.gt.this%endphen)exit
  k=k+1
  this%pheno(:,k)=rtmp
  this%presentpheno(k)=i8tmp
  this%effect(:,k)=itmp
 enddo
 close(un)
 this%effect=this%effect-1

 !read residual variance file
 allocate(this%resmati(this%ntrait,this%ntrait))
 open(newunit=un,file=this%resfile,status='old',action='read')
 do i=1,this%ntrait
  read(un,*)rtmp(1:i)
  this%resmati(i,1:i)=rtmp(1:i)
  this%resmati(1:i,i)=rtmp(1:i)
 enddo
 close(un)
 

 !Max combi
 maxcombi=maxval(this%presentpheno)
 allocate(this%tres(maxcombi))
 do i8tmp=0,maxcombi
  if(any(this%presentpheno.eq.i8tmp))then
   allocate(this%tres(i8tmp)%mat(this%ntrait,this%ntrait))
   this%tres(i8tmp)%mat=this%resmati
   do i=1,this%ntrait
    if(ibits(i8tmp,i-1,1).eq.0)then
     this%tres(i8tmp)%mat(:,i)=0._real8
     this%tres(i8tmp)%mat(i,:)=0._real8
     this%tres(i8tmp)%mat(i,i)=1._real8
    endif
   enddo
   call invsym(this%tres(i8tmp)%mat,this%ntrait,this%ntrait)
   do i=1,this%ntrait
    if(ibits(i8tmp,i-1,1).eq.0)this%tres(i8tmp)%mat(i,i)=0._real8
   enddo
  endif
 enddo

! do i8tmp=0,maxcombi
!  if(allocated(this%tres(i8tmp)%mat))write(this%unlog,*)'ress ',this%tres(i8tmp)%mat
! enddo

 !invert resmati
! call invsym(this%resmati,this%ntrait,this%ntrait)
! write(this%unlog,*)'bbbb ',this%resmati

end subroutine

!**MULTIPLICATION BY V
#if (COARRAY==1)
subroutine multbyv_mvlr(this,x,y,starteq,endeq)
 !y=M*x=X'*Ri*X*x
 class(mvlr),intent(inout)::this
 integer(kind=int4),intent(in)::starteq,endeq
 real(kind=real8),intent(in)::x(:)
 real(kind=real8),intent(inout)::y(:)

 integer(kind=int4)::i,j,k,address
 real(kind=real8),allocatable::t1(:),t2(:),xmat(:,:),ri(:,:)
 real(kind=real8),allocatable,save::ylong(:)[:]
 !$ real(kind=real8)::time1,time2
 
 !$ time1=omp_get_wtime()
 if(.not.allocated(ylong))then
  !$ time2=omp_get_wtime()
  allocate(ylong(this%neq)[*])
  write(this%unlog,'(a)')' allocation of ylong'
  !$ write(this%unlog,'(a,f0.3)')'  Elapsed time(s) alloc: ',omp_get_wtime()-time2
 endif


 ylong=0_real8
 
 sync all

 y=0_real8

 allocate(t1(this%ntrait),t2(this%ntrait))
 allocate(xmat(this%ntrait,this%neffect),ri(this%ntrait,this%ntrait))

 !y=X'*Ri*X*x
 do i=1,this%nsubphen
  !get X
  call getX(this,xmat,i)
  !get Ri
  ri=this%tres(this%presentpheno(i))%mat
  !1) t1=X*x
  t1=0.d0
  do j=1,this%neffect
   address=this%effect(j,i)
   do k=1,this%ntrait
    t1(k)=t1(k)+xmat(k,j)*x(address+k)
   enddo
  enddo
  !2) t2=Ri*t1
  t2=0.d0
  do j=1,this%ntrait
   do k=1,this%ntrait
    t2(j)=t2(j)+ri(k,j)*t1(k)
   enddo
  enddo
  !3) y=X'*t2
  do j=1,this%neffect
   address=this%effect(j,i)
   do k=1,this%ntrait
    ylong(address+k)=ylong(address+k)+xmat(k,j)*t2(k)
   enddo
  enddo
 enddo

 deallocate(t1,t2,xmat,ri)
 
 sync all

 !update y
 do i=1,num_images()
  y=y+ylong(starteq:endeq)[i]
 enddo

 sync all
 
 !$ write(this%unlog,'(a,f0.3)')'  Elapsed time(s) multbyv: ',omp_get_wtime()-time1

end subroutine

#else
!**PRECONDITIONER
function computediag_vect_mvlr(this) result(array)
 class(mvlr),intent(inout)::this
 real(kind=real8),allocatable::array(:)

 integer(kind=int4)::i,j,k,l,m,address1,address2
 real(kind=real8)::val
 real(kind=real8),allocatable::xmat(:,:),ri(:,:)
 type(coosparse)::sparse

 allocate(xmat(this%ntrait,this%neffect),ri(this%ntrait,this%ntrait))

 if(allocated(array))deallocate(array)
 
 sparse=coosparse(this%neq,this%neq,int(this%neq,8),lupper=.true.)

 do i=1,this%nsubphen
  !get X
  call getX(this,xmat,i)
  !get Ri
  ri=this%tres(this%presentpheno(i))%mat
  
  do j=1,this%neffect
   address1=this%effect(j,i)
!   do k=1,this%neffect
    k=j
    address2=this%effect(k,i)
    do l=1,this%ntrait
!     do m=1,this%ntrait
      m=l
      val=xmat(l,j)*ri(m,l)*xmat(m,k)
      call sparse%add(address1+l,address2+m,val)
!     enddo
    enddo
!   enddo
  enddo
 enddo
 
 call sparse%printstats() 

 array=sparse%diag()

end function

#endif

!PRIVATE
subroutine getX(this,xmat,iphen)
 type(mvlr),intent(inout)::this
 integer(kind=int4),intent(in)::iphen
 real(kind=real8),intent(inout)::xmat(:,:)

 integer(kind=int4)::i,j

 xmat=1_real8
 do j=1,this%neffect
  do i=1,this%ntrait
   if(.not.this%model(i,j))xmat(i,j)=0_real8
  enddo
 enddo

end subroutine

subroutine invsym(a,n,m)
 !Inverse of a real symmetric positive definite matrix of size n*n declared as m*m 
 integer(kind=int4),intent(in)::n,m
 real(kind=real8),intent(inout)::a(:,:)

 integer(kind=int4)::info,i,j
 
 call dpotrf('L',n,a,m,info)
 if (info.ne.0) then
  print*,'ERROR dpotrf info: ',info
  stop
 endif
 call dpotri('L',n,a,m,info)
 if (info.eq.0) then
 !$omp parallel 
 !$omp do schedule(auto)
  do i=1,n
   do j=i+1,n
    a(i,j)=a(j,i)
   enddo
  enddo
 !$omp enddo
 !$omp end parallel 
 else
  print*,'ERROR dportri info: ',info
  stop
 endif

end subroutine

end module
