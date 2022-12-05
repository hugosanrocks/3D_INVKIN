!copyright(c) 2016 seiscope consortium
!Author: Pengliang Yang, Univ. Grenoble Alpes
!E-mail: pengliang.yang@univ-grenoble-alpes.fr
!Adapted version from Jon Claerbout
!===========================================================================
program smooth2d
  implicit none

  integer::n1,n2,r1,r2, repeat
  real,dimension(:,:),allocatable::mod0,mod1
  integer::i1,i2,i
  character(len=80)::inputfile,outputfile
  real::model(41,648)
  integer::ii,j,k,l,iunit1,iunit2


  allocate(mod0(18,36))
  allocate(mod1(18,36))

  iunit1=22
  iunit2=23
  open(iunit1,file='prior_model.dat',status='old')
  open(iunit2,file='prior_model_smooth.dat',status='unknown')
  do ii=1,41
   read(iunit1,*) model(ii,:)
   l=1
   do j=1,18
    do k=1,36
     mod0(j,k) = model(ii,l)
     l=l+1
    enddo
   enddo
  

  !read(*,*) inputfile,outputfile
  !read(*,*) n1,n2,r1,r2,repeat !size of data, smoothing radius,repeating times
!  print *,inputfile,outputfile,n1,n2,r1,r2,repeat
  n1=18
  n2=36
  r1=1
  r2=2
  repeat=2

  !allocate(mod0(n1,n2))
  !allocate(mod1(n1,n2))

  !input to  be smoothed
  !open(10,file=inputfile,access='direct',recl=4*n1*n2)
  !read(10,rec=1) mod0
  !close(10)
  do i=1,repeat
     do i2=1,n2
        call triangle(r1,n1,mod0(:,i2),mod1(:,i2))
     enddo
     do i1=1,n1
        call triangle(r2,n2,mod1(i1,:),mod0(i1,:))
     enddo
  enddo
  l=1
  do j=1,18
   do k=1,36
     if (j .gt. 15) then
      model(ii,l) = model(ii,l)
     else
      model(ii,l) = mod0(j,k)
     endif
     l=l+1
   enddo
  enddo
  write(iunit2,*) model(ii,:)
  enddo
  close(iunit2)
  close(iunit1)
  !smoothed version of the input
  open(10,file=outputfile,access='direct',recl=4*n1*n2,status='replace')
  write(10,rec=1) mod0
  close(10)


  deallocate(mod0)
  deallocate(mod1)

  call exit(0)
end program smooth2d

!==================================================================
subroutine boxconv(nbox, nx, xx, yy)
  implicit none
  integer,intent(in):: nx,nbox
  integer::i,ny
  real,dimension(nx),intent(in)::xx
  real,dimension(nx+nbox-1),intent(out)::yy
  real,dimension(:),allocatable::bb

  allocate(bb(nx+nbox))
  if (nbox < 1 .or. nbox > nx) then
     write(0,*) "boxconv: error in the length of input!"
  endif
  ny = nx+nbox-1
  do i= 1, ny
     bb(i) = 0.
  enddo
  bb(1) = xx(1)
  do i= 2, nx
     bb(i) = bb(i-1) + xx(i)  ! make B(Z) = X(Z)/(1-Z)
  enddo
  do i= nx+1, ny
     bb(i) = bb(i-1)
  enddo
  do i= 1, nbox
     yy(i) = bb(i)
  enddo
  do i= nbox+1, ny
     yy(i) = bb(i) - bb(i-nbox) ! make Y(Z) = B(Z)*(1-Z**nbox)
  enddo
  do i= 1, ny
     yy(i) = yy(i) / nbox
  enddo
  deallocate(bb)
end subroutine boxconv

!=================================================================
!apply triangle filter to input xx: yy=triangle*xx
subroutine triangle(nbox,nd,xx,yy)
  implicit none

  integer, intent(in)::nbox,nd
  integer::i,np,nq
  real,dimension(nd),intent(in)::xx
  real,dimension(nd),intent(out)::yy
  real,dimension(:),allocatable::pp,qq

  allocate(pp(nd+nbox-1))
  allocate(qq(nd+2*nbox-2))
  call boxconv(nbox,nd,xx,pp)
  np=nbox+nd-1
  call boxconv(nbox,np,pp,qq)
  nq=nbox+np-1
  do i=1,nd
     yy(i)=qq(i+nbox-1)
  enddo
  do i=1,nbox-1
     yy(i)=yy(i)+qq(nbox-i) !fold back near end
  enddo
  do i=1,nbox-1
     yy(nd-i+1)=yy(nd-i+1)+qq(nd+(nbox-1)+i) !fold back far end
  enddo
  deallocate(pp)
  deallocate(qq)
end subroutine triangle
