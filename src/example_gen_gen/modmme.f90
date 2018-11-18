module modmme
 use iso_fortran_env
 !$ use omp_lib
 use modsparse
 use modcoeff,only:gen_coeff
 use modmvlr

 implicit none
 private
 public::mme

#if (COARRAY==1)
 type,extends(gen_coeff)::mme
#else
 type::mme
#endif
  private
  integer(kind=int32)::unlog=6
  integer(kind=int32)::neq
  integer(kind=int32)::nsparse
  character(len=30),allocatable::namesparsefile(:)
  type(mvlr)::lr
  type(crssparse),allocatable::mat(:)
 contains
  private
  procedure,public::init=>init_mme
  procedure,public::getneq
  procedure,public::multbyv=>multbyv_mme
 end type

contains

!PUBLIC
!**INIT
subroutine init_mme(this,paramfile,optparfile,startrow,endrow,unlog)
 class(mme),intent(inout)::this
 integer(kind=int32),intent(in),optional::startrow,endrow
 integer(kind=int32),intent(in),optional::unlog
 character(len=*),intent(in)::paramfile
 character(len=*),intent(in),optional::optparfile

 integer(kind=int32)::i,un,io
 type(crssparse)::mat

 this%unlog=6
 if(present(unlog))this%unlog=unlog

 call this%lr%init(paramfile,this%unlog)
 
 this%neq=this%lr%getneq()

 if(present(optparfile))then
  if(.not.present(startrow).or..not.present(endrow))then
   write(this%unlog,'(a)')' ERROR: startrow and endrow must be present with optparfile!'
   error stop
  endif
  write(this%unlog,'(/2a/)')' Optional parameter file: ',trim(optparfile)
  open(newunit=un,file=optparfile,status='old',action='read')
  read(un,*,iostat=io)this%nsparse
  if(io.ne.0)error stop
  if(this%nsparse.gt.0)then
   allocate(this%namesparsefile(this%nsparse))
   this%namesparsefile=''
   allocate(this%mat(this%nsparse))
   do i=1,this%nsparse
    read(un,*,iostat=io)this%namesparsefile(i)
    if(io.ne.0)exit
    mat=crssparse(this%namesparsefile(i))
    if(this%neq.ne.mat%getdim(1).or..not.mat%lsquare())then
     write(this%unlog,'(2a)')' ERROR: invalid sparse matrix: ',this%namesparsefile(i)
     error stop
    endif
    this%mat(i)=mat%submatrix(startrow,endrow,1,this%neq)
    call this%mat(i)%setoutputunit(this%unlog)
    call this%mat(i)%printstats()
    call mat%destroy()
 
    !for testing
    !mat=this%mat(i)%submatrix(1,10,1,10)
    !call mat%setoutputunit(this%unlog)
    !call mat%printsquare()
    !call mat%destroy()

   enddo
  endif
  close(un)
 endif

end subroutine

!**GET #EQUATIONS
function getneq(this) result(neq)
 class(mme),intent(in)::this
 integer(kind=int32)::neq

 neq=this%neq
end function

!**MULTIPLICATION (FUNCTION)
subroutine multbyv_mme(this,x,y,starteq,endeq)
 class(mme),intent(inout)::this
 integer(kind=int32),intent(in)::starteq,endeq
 real(kind=real64),intent(in)::x(:)
 real(kind=real64),intent(inout)::y(:)

 integer(kind=int32)::i

 call this%lr%multbyv(x,y,starteq,endeq)  

 sync all

 do i=1,this%nsparse
  call this%mat(i)%multbyv(1._real64,'n',x,1._real64,y)
 enddo

 sync all

end subroutine

end module
