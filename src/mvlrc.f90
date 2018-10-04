program mvlrc
 use modkind
 use modmvlr
 implicit none
 integer(kind=int4)::neq
 integer(kind=int4)::nsize,starteq,endeq
 real(kind=real8),allocatable::x(:)
 real(kind=real8),allocatable::y(:)
 type(mvlr)::reg

 call reg%init('parammvlr.dat',500+this_image())

 neq=reg%getneq()

 allocate(x(neq))
 x=1_real8

 
 nsize=int(real(neq)/num_images())
 starteq=nsize*(this_image()-1)+1
 endeq=starteq+nsize-1
 if(this_image().eq.num_images())endeq=neq


 allocate(y(endeq-starteq+1))

 !y=X'*Ri*X*x
 call reg%multbyv(x,y,starteq,endeq)

 write(500+this_image(),*)'a ',starteq,' b ',endeq,' aaaa ',y

 x=0.5
 call reg%multbyv(x,y,starteq,endeq)
 write(500+this_image(),*)'a ',starteq,' b ',endeq,' aaaa ',y

contains


end program
