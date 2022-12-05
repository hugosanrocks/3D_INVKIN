       program test

       implicit none
  include 'mpif.h'

  integer, parameter :: ncell_min= 40
  integer, parameter :: ncell_max= 700

  real, parameter :: mean = 2.8 ! Mean of the Uniform prior on response value 
  real, parameter :: theta = 0.6

  real t1,t2
  integer, parameter :: nt_max=20000
  integer, parameter :: nv_max=10000
  integer, parameter :: nmax=3*nt_max+ncell_max


  !Bounds of the 2D region
  real, parameter :: longmin = 0., longmax = 34.5
  real, parameter :: latmin = 0., latmax = 16.5

  double precision  :: eps = 0.00000001
  integer           nnn(ncell_max+1)
  integer           nnlist(nmax)
  integer           ntwork(nmax)
  integer           vertices(3,nt_max)
  integer           neighbour(3,nt_max)
  integer           worki1(nv_max)
  integer           worki2(nv_max)
  integer           worki3(nv_max)
  integer           nt

  real , EXTERNAL    ::    gasdev,logprior,ran3
  integer ncell
  double precision Voro(3,ncell_max)

  logical	    ldummy(ncell_max)
  integer i, j, p

  !for MPI
  integer ra,rank, nbproc, ierror, tag, status(MPI_STATUS_SIZE)

  CALL cpu_time(t1)! Tic. start counting time 

  rank = 0
  ra=rank

ncell=int(ncell_min+ran3(ra)*(ncell_max-ncell_min))
ncell=500

p=1
!* We place them randomly
DO i=1,17!ncell
 ! Initial location of the cells (V(:,1:2))
!   Voro(1:2,i) = [longmin+ran3(ra)*(longmax-longmin),latmin+ran3(ra)*(latmax-latmin)]
do j=1,35
 if ( p .eq. 1 ) then
  voro(1:2,i) = [ real(j)-0.75,real(i)-0.5 ]
 else
  voro(1:2,i) = [ real(j)-0.25,real(i)-0.5 ]
 endif
write(77,*) Voro(1:2,i)
enddo
 if (p .eq. 1) then
  p=2
 else
  p=1
 endif
ENDDO

! compute the Voronoi grid
call delaun (Voro(1:2,:),ncell,neighbour,vertices,nt,nt_max,&
                      worki1,worki2,worki3,eps,nv_max,&
                      0,ldummy,0,0,0)

print *, nt
do i=1,nt
 write(78,*) vertices(1:3,i)
enddo

!call build_nv(ncell,vertices,nt,ncell_max,nmax,&
!              neighbour,nnn,nnlist,ntwork)

!call MPI_FINALIZE(ierror)! Terminate the parralelization


CALL cpu_time(t2) !toc stop counting time and display
if (rank==0) write(*,*)'time taken by the code was',t2-t1,'seconds'


       end

!-------------------------------------------------------------------
!                                               
!       Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
!      write(*,*)' idum ',idum
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

