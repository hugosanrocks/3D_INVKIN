PROGRAM yoffee_delay

  REAL, DIMENSION(:,:), ALLOCATABLE :: STF
  REAL t, dt, dist, tac, stk_step, dip_step
  REAL rt, tfin, tshift, tac1
  INTEGER i, j, nstk, ndip, tsam


  tfin = 2.
  tsam = 2. / 0.04
  nstk = 50
  ndip = 20
  stk_step = 0.4
  rt   = 1.
  vr   = 2.9


  allocate(STF(tsam,nstk))


  tac1 = (nstk_step / 2) / vr
  dt = 0.04 !tfin / real(tsam,4)

  t=0.
  do i=1,tsam
  STF(:,1) = yoffe(t,dt,2.,0.2)
  t = t + dt
  print *, STF(i,1)
  enddo


  tac = dt
 
  dist = 0.0

  do i = 1,nstk-1

    tshift = dist / vr
    !print *, tshift

    dist = dist + stk_step
  enddo


  deallocate(STF)

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
  REAL, PARAMETER :: Dm = 2.0      	! Total slip
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
