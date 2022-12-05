  PROGRAM yoffee_delay
  
  IMPLICIT NONE
  REAL, DIMENSION(:), ALLOCATABLE :: dist, tshift, tssam, tinter, stsig
  REAL, DIMENSION(:,:), ALLOCATABLE :: STF, STFI, STFT
  REAL t, dt, stk_step, dip_step
  REAL rt, tfin, tac1, vr, re, me
  INTEGER, DIMENSION(:), ALLOCATABLE :: tsamint, tsamfin
  INTEGER i, j, k, m, nstk, ndip, tsam, sigsam

  INTEGER j8, j9
  REAL const, slope, corr, x(2), y(2), intval

  REAL yoffe

  tfin = 25.
  tsam = 2048
  nstk = 50
  ndip = 20
  stk_step = 0.4
  rt   = 1.0
  vr   = 2.9
  dt = 0.01227
  tac1 = 0.2

  sigsam = int(((rt / dt)+1),4)


  allocate(STF(tsam,nstk),tsamfin(nstk),stsig(sigsam),STFI(tsam,nstk))
  allocate(dist(nstk),tshift(nstk),tssam(nstk),tinter(nstk),tsamint(nstk))
  allocate(STFT(tsam,nstk*ndip))


  dist(1) = 0.
  do i = 2,nstk
     dist(i) = dist(i-1) + stk_step
  enddo
  tshift(:) = dist(:) / vr
  tssam(:) = tshift(:) / dt

  do i = 1,nstk
   tsamint(i) = int(tssam(i),4)
   tinter(i) = tssam(i) - real(int(tssam(i),4))
   tsamint(i) = tsamint(i) + 1
   tsamfin(i) = tsamint(i) + sigsam
  enddo


  STF(:,:) = 0.

  t = 0
  do i = 1,tsam - 1
  STF(i,1) = yoffe(t,dt,rt,tac1)
  t = t + dt
  enddo

  stsig(:) = STF(1:sigsam,1)
  STFI(:,1) = STF(:,1)

  do i= 2,nstk
   m = tsamint(i)
   k = tsamfin(i)
   STF(m:k,i) = stsig(:)
  enddo

!  do i = 1,tsam
!   write(33,*) STF(i,:)
!  enddo

  J8 = 2
  J9 = 0

  t = 0
  do j = 2,nstk
  do i = 1,tsam-1
   y(1) = STF(i,j)
   y(2) = STF(i+1,j)
   x(1) = t - ((1.-tinter(j))*dt)
   x(2) = x(1) + dt

   if ((nstk .eq. 5) .or. (j .eq. 9)  ) then

   STFI(i,j) = STF(i,j)

   else

   call linreg(J8,J9,CONST,SLOPE,CORR,X,Y)
   intval = (t*slope) + const
   STFI(i,j) = intval

   endif

   t = t + dt
  enddo
  enddo

  !Need to check problem at this times
  STFI(:,5) = STF(:,5)
  STFI(:,6) = STF(:,6)
  STFI(:,10) = STF(:,10)


  do i=1,nstk
   do j=1,ndip
    m = i+((j-1)*nstk)
    STFT(:,m) = STFI(:,i)
   enddo
  enddo


  OPEN(33,FILE='yoffe.src',status='unknown')
  do i = 1,tsam
   write(33,*) STFT(i,:)
  enddo
  close(33)







  deallocate(STFT)
  deallocate(dist,tshift,tssam,tinter,tsamint,STF,stsig,STFI)
  ENDPROGRAM yoffee_delay

! ++++++++++++++++++++++++++++++
!
  FUNCTION yoffe(t,dt,rre,tacc)
!
! Translation to Fortran90 of the Matlab routine provided by 
! Elisa Tinti to generate the regularized Yoffe slip-rate 
! function (see Tinti et al., 2005, BSSA) 
!
! By Victor M. Cruz-Atienza (2007)
!
  IMPLICIT NONE
!
! Input Parameters
!
  REAL, PARAMETER :: Dm = 1.0      	! Total slip
!                    dt = 0.005, &	! Time step
!                    rre = 2.0, &   	! Effective rise-time
!                    tacc = 0.2     	! Time AOF acceleration
!
! Variable Declaration
!
  INTEGER :: n,nt,i
  REAL :: dt,rre,tacc
  REAL :: ss,rr,t
  REAL :: CO,CC2,DD2
  REAL :: pi
  REAL :: AA,BB,CC,DD
  REAL :: AA1,BB1,CC1,DD1
  REAL :: AA3,BB3,CC3,DD3
  REAL :: AA4,BB4
  REAL :: yoffe
!
  ss=tacc/1.27   	! tau_s = half duration of triangular function
  rr=rre-2*ss    	! rr = rise time of singular yoffe
! 
! Effective Rise Time rre
! rre=rr+2*ss
!
  pi=4.*ATAN(1.)
  CO=2./(pi*ss**2.*rr)*Dm
  nt=NINT((rr+2.*ss)/dt)+1
!
  IF (t <= rr+2*ss) THEN
    IF (t <= ss) THEN
!
      CC2=sqrt(t*(rr-t))*(1./2.*t+1./4.*rr)+asin(sqrt(t/rr))* &
              (t*rr-rr**2)-3./4.*rr**2*atan(sqrt((rr-t)/t))
      DD2=3./8.*pi*rr**2
      yoffe=(CC2+DD2)
!
    ELSE IF (t > ss .AND. t <= rr+ss) THEN
!
      IF (rr > 2*ss) THEN
        IF (t <= 2*ss) THEN
!
          AA1=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
              (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
               3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB1=-3./8.*pi*rr**2
          CC1=sqrt(t*(rr-t))*(1./2.*t+1./4.*rr)+(t*rr-rr**2)* &
              asin(sqrt(t/rr))-3./4.*rr**2*atan(sqrt((rr-t)/t))
          DD1=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
              (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+3./4.*rr**2* &
              atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe=(AA1+BB1+CC1+DD1)
!
        ELSE IF (t > 2*ss .AND. t <= rr) THEN
!
          AA=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
             (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
             3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB=sqrt((t-2*ss)*(rr-t+2*ss))*(-ss+1./4.*rr+1./2.*t)+ &
             (-2*ss*rr+t*rr-rr**2)*asin(sqrt((t-2*ss)/rr))- &
             3./4.*rr**2*atan(sqrt((rr-t+2*ss)/(t-2*ss)))
          CC=sqrt(t*(rr-t))*(1./2.*t+1./4.*rr)+(t*rr-rr**2)* &
             asin(sqrt(t/rr))-3./4.*rr**2*atan(sqrt((rr-t)/t))
          DD=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
             (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+3./4.*rr**2* &
             atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe=(AA+BB+CC+DD)
!
        ELSE IF (t > rr .AND. t <= rr+ss) THEN
!
          AA3=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
              (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB3=sqrt((t-2*ss)*(rr-t+2*ss))*(-ss+1./4.*rr+1./2.*t)+ &
              (-2*ss*rr+t*rr-rr**2)*asin(sqrt((t-2*ss)/rr))- &
              3./4.*rr**2*atan(sqrt((rr-t+2*ss)/(t-2*ss)))
          CC3=pi*0.5*(t*rr-rr**2)
          DD3=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
              (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe=(AA3+BB3+CC3+DD3)
!
        END IF
!
      ELSE IF (rr <= 2*ss) THEN
!
        IF (t <= rr) THEN
!
          AA1=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
              (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB1=-3./8.*pi*rr**2
          CC1=sqrt(t*(rr-t))*(1./2.*t+1./4.*rr)+(t*rr-rr**2)* &
              asin(sqrt(t/rr))-3./4.*rr**2*atan(sqrt((rr-t)/t))
          DD1=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
              (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe= (AA1+BB1+CC1+DD1)
!
        ELSE IF (t <= 2*ss .AND. t > rr) THEN
!
          AA=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
             (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
             3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB=-3./8.*pi*rr**2
          CC=pi*0.5*(t*rr-rr**2)
          DD=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
             (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
             3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe=(AA+BB+CC+DD)
!
        ELSE IF (t <= (rr+ss) .AND. t > 2*ss) THEN
!
          AA3=sqrt((t-ss)*(rr-t+ss))*(3./2.*ss-1./4.*rr-1./2.*t)+ &
              (2*ss*rr-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          BB3=sqrt((t-2*ss)*(rr-t+2*ss))*(-ss+1./4.*rr+1./2.*t)+ &
              (-2*ss*rr+t*rr-rr**2)*asin(sqrt((t-2*ss)/rr))- &
              3./4.*rr**2*atan(sqrt((rr-t+2*ss)/(t-2*ss)))
          CC3=pi*0.5*(t*rr-rr**2)
          DD3=sqrt((t-ss)*(rr-t+ss))*(-1./2.*t-1./2.*ss-1./4.*rr)+ &
              (-t*rr+rr**2)*asin(sqrt((t-ss)/rr))+ &
              3./4.*rr**2*atan(sqrt((rr-t+ss)/(t-ss)))
          yoffe=(AA3+BB3+CC3+DD3)
!
        END IF
      END IF
!
    ELSE IF (t > rr+ss) THEN
!
      AA4=pi*0.5*(rr*(2*ss-t+rr))
      BB4=sqrt((t-2*ss)*(rr-t+2*ss))*(-ss+1./4.*rr+1./2.*t)+ &
          (-2*ss*rr+t*rr-rr**2)*asin(sqrt((t-2*ss)/rr))- &
          3./4.*rr**2*atan(sqrt((rr-t+2*ss)/(t-2*ss)))
      yoffe=(AA4+BB4)
!
    END IF
    yoffe=yoffe*CO
!
  ELSE
    yoffe=0.0
  END IF
!
  END FUNCTION yoffe
!
! ++++++++++++++++++++++++++++++
!
