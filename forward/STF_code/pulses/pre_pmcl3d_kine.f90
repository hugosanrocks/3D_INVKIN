!

PROGRAM pre_pmcl3d_kine
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: ifpoint=20000000
! REAL, PARAMETER :: mu=2.83E10
! REAL, PARAMETER :: mu=3.2E10    ! Vp=6.0, Vs=3.464, rho=2.67
! REAL, PARAMETER :: mu=3.328E10  ! Vp=6682.0, Vs=3858.0, rho=2908.0
!  REAL, PARAMETER :: mu=7.681E10  
  INTEGER, PARAMETER :: nswid=381, nslen=251  ! Dimensions of slip input file
!
! Variable Declaration
!
  INTEGER :: i,j,it,ifa
  INTEGER :: ntall,nsubf,nini
  INTEGER :: kstf,istf,iread
  INTEGER :: inp3d,icont,itmp,icont1
  INTEGER :: ixn,iyn,izn
  INTEGER :: nlength,nwide
  INTEGER :: fflag,sflag,iout,dipole
  REAL :: mu
  REAL :: dx,dt,pi,mw
  REAL :: tr,rvel,tt,t0
  REAL :: rakes,raked,str,dip
  REAL :: axx_s,ayy_s,azz_s,axz_s,ayz_s,axy_s
  REAL :: axx_d,ayy_d,azz_d,axz_d,ayz_d,axy_d
  REAL :: xx_s,yy_s,zz_s,xz_s,yz_s,xy_s
  REAL :: xx_d,yy_d,zz_d,xz_d,yz_d,xy_d
  REAL :: area,rtime,dist
  REAL :: cnx,cny,cnz
  REAL :: cx,cy,cz
  REAL :: areatot,m0,avgslip,srarea,areamuctte
  REAL :: lof,hif
  REAL :: rre,tacc
  REAL :: trash, vp, vs, rho
  REAL :: velo(350,70)
  REAL :: sltmp
  CHARACTER(LEN=30) :: file1,file2,file3,ftype
!
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: coord3d,nucle,coord2d
  REAL, DIMENSION(:,:), ALLOCATABLE :: vnormal,fault_mech, &
                                       rsrs,rsrd
! REAL, DIMENSION(:), ALLOCATABLE :: srs,srd,srtot,m0t,fslip,areamu
  REAL, DIMENSION(:), ALLOCATABLE :: srs,srd,srtot,m0t,areamu,fslip
! REAL, DIMENSION(:) :: fslip(95631)
! REAL, DIMENSION(:) :: fslip(100000)
!
  PRINT *
  WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*,*) '               PROGRAM: PRE_PMCL3D_KINE'
  WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT *
!
! Input Files
!
  OPEN(32,FILE='pre_pmcl3d_kine.in',STATUS='old')
  OPEN(33,FILE='fault_geometry.in',STATUS='old')
  OPEN(90,FILE='ModelSvelo_nvm.dat', STATUS='old')
!
! Read Input Data
!
  DO i=1,4
    READ(32,*)
  END DO
  READ(32,*) iread
  READ(32,*) file1,file2
  READ(32,*) sflag
  READ(32,*) file3
  READ(32,*) ntall
  READ(32,*) dx,dt
  READ(32,*) fflag,ftype
  READ(32,*) lof,hif
  READ(32,*) rvel
  READ(32,*) mw
  READ(32,*) iout
  READ(32,*) dipole
  READ(32,*) kstf
  DO i=1,5
    READ(32,*)
  END DO
  READ(32,*) istf
  READ(32,*) tr
  DO i=1,5
    READ(32,*)
  END DO
  READ(32,*) rre
  READ(32,*) tacc
!
! Real Velocity Model
!
DO i=1,350
  DO j=1,70
    READ(90,*) trash, trash, velo(i,j)
  END DO
END DO
!
  pi=4.0*ATAN(1.0)
! area=dx**2.0
  area=3e3**2.0
!
! Precomputed Slip-Rate Functions Input Files
!
  IF (iread == 1) THEN
!   OPEN(17,FILE=file1)
    OPEN(18,FILE=file2)
  END IF
!
! Precomputed Final Slip Functions Input File
!
! (The one column final slip ascii file must be written first along
!  the fault-dip index and then along the fault-strike index)
!
  IF (sflag == 1.AND.iread == 0) THEN
    OPEN(19,FILE=file3)
  END IF
!
! Output Files
!
  OPEN(31,FILE='rtimes.dat',STATUS='unknown')
  OPEN(34,FILE='coord2D.dat',STATUS='unknown')
  OPEN(35,FILE='coord3D.dat',STATUS='unknown')
  OPEN(36,FILE='fault_mech.dat',STATUS='unknown')
  IF (iout == 1) THEN
    OPEN(37,FILE='momrate',STATUS='unknown')
  END IF
  OPEN(38,FILE='slip-rate',STATUS='unknown')
  OPEN(39,FILE='slip',STATUS='unknown')
!
! Array Memory Allocation
!
  ALLOCATE ( coord3d(ifpoint,3),nucle(ifpoint,4), &
             vnormal(ifpoint,3), fault_mech(ifpoint,3), &
             srs(ntall), srd(ntall), srtot(ntall), &
             m0t(ntall), coord2d(ifpoint,2), areamu(ifpoint), &
             fslip(ifpoint) )
!
! Create Fault Geometry
!
  WRITE(*,*) 'Computing Fault Geometry:'
  CALL geometry(ifpoint,nsubf,inp3d,dx,coord3d, &
                vnormal,nucle,fault_mech,ixn,iyn,izn, &
                nlength,nwide,coord2d)
!
! Read Velocity at the fault 
!
  DO j=1,nlength
    vs=velo(coord2d(j,1), coord2d(j,2))*1000
    vp=sqrt(3.00)*vs
    rho=(0.32*vp/1000+0.77)*1000
    areamu(j)=rho*vp**2*area
  END DO
!
! Read or Compute Final Slip 
!
  !write(*,*) nswid,nslen,nwide,nlength
  IF (iread == 0) THEN
    IF (sflag == 1) THEN
      !DO i=1,nswid*nslen
      DO i=1,nwide*nlength
        READ(19,*) fslip(i)
      END DO
    ELSE
      areatot=nlength*nwide*area
      m0=10.0**(1.5*mw+16.1)*1E-5
      avgslip=m0/(mu*areatot)
      fslip(:)=avgslip
      WRITE(*,FMT=51) mw,avgslip
      IF (avgslip > 15.0) WRITE(*,*) 'WARNING!!! You should decrease Mw'
      PRINT *
    END IF
  END IF

! fslip(:)=fslip(:)*40.0    ! To scales SSE slip (remove)

!
  51 FORMAT(' The constant final slip for the Mw=',F3.1 &
           ,' event is:',F5.1,' m')
!
! Read or Compute Slip-Rate Funtions
!
  srtot(:)=0.0
!
  IF (iread == 1) THEN
    ALLOCATE ( rsrs(nsubf,ntall), rsrd(nsubf,ntall) )
    DO ifa=1,nsubf
      DO it=1,ntall
!       READ(17,*) rsrs(ifa,it)
        READ(18,*) rsrd(ifa,it)
      END DO
    END DO
!   rsrs(:,:)=0.0      ! Remove this line !!!!!!
  ELSE
    tt=0.0
    DO it=1,ntall
      IF (kstf == 1.AND.istf == 1) srtot(it)=trian(tr,tt,t0,dt,pi)
      IF (kstf == 1.AND.istf == 2) srtot(it)=gauss(tr,tt,t0,dt,pi)
      IF (kstf == 2) srtot(it)=yoffe(tt,dt,rre,tacc)
      tt=tt+dt
    END DO
  END IF
!
! Align Computed Slip-Rate Function with t0
!
  itmp=0
  IF (itmp == 1) THEN
  icont=2
  DO it=2,ntall
    IF (srtot(it) > 1E-15) THEN
      srtot(icont)=srtot(it)
      IF (srtot(icont) >= 0.01.AND.nini == 0) nini=icont
      icont=icont+1
    END IF
  END DO
  DO it=icont,ntall
    srtot(it)=0.0
  END DO
  END IF
!
! Filter Slip-Rate Function (If asked for)
!
  IF (fflag == 1) CALL xapiir(srtot,ntall,'BU',0.0,0.0,4, &
                              ftype,lof,hif,dt,1)
!
! Normalize Slip-Rate Function (Its integral equals one)
!
  srarea=0.0
  DO it=1,ntall
    srarea=srarea+srtot(it)*dt
  END DO
  srtot(:)=srtot(:)/srarea
!
! Compute Slip Function (for visualization)
!
  tt=0.0
  m0t(:)=0.0
!
  WRITE(38,FMT=42) tt,srtot(1) 
  WRITE(39,FMT=42) tt,m0t(1) 
  DO it=2,ntall
    tt=tt+dt
    m0t(it)=m0t(it-1)+dt*srtot(it)
    WRITE(38,FMT=42) tt,srtot(it) 
    WRITE(39,FMT=42) tt,m0t(it) 
  END DO
!
! Convert Fault Mechanism Angles (strike, dip, rake) into Radians
!
  fault_mech(:,:)=fault_mech(:,:)*pi/180.0
!
! Coordinates Nucleation Point
!
  cnx=(ixn-1)*dx
  cny=(iyn-1)*dx
  cnz=(izn-1)*dx
!
! ++++++++++++++++++++++++++++++++++++++
!
  WRITE(*,*) 'Computing Fault Moment-Rate Tensor History:'
  PRINT *
!
! Loop Over Subfaults
!
  icont1=0
!
  DO ifa=1,nsubf
!
    IF (MOD(ifa,nwide) == 0 .OR. ifa == 1) THEN
      icont1=icont1+1
    END IF
!
    IF (MOD(ifa,500) == 0.AND.iout == 1) WRITE(*,FMT=50) ifa,nsubf
    50 FORMAT('      Subfault number',I6,' of',I6)
!
! Write Current Subfault Indices
!
    IF (iout == 1) THEN
      WRITE(37,FMT=40) coord3d(ifa,1),coord3d(ifa,2),coord3d(ifa,3)
    END IF
!
! Decompose Slip-Rate Functions into Along-Strike and 
! Along-Dip Components and Scale to Match the Final Slip
!
    IF (iread == 1) THEN
!     srs(:)=rsrs(ifa,:)
      srd(:)=rsrd(ifa,:)
    ELSE IF (iread == 0.AND.dipole == 0) THEN
      !srs(:)=srtot(:)*ABS(COS(fault_mech(ifa,3)))*fslip(ifa)*0.5
      !srd(:)=srtot(:)*ABS(SIN(fault_mech(ifa,3)))*fslip(ifa)*0.5
      srs(:)=srtot(:)*ABS(COS(fault_mech(ifa,3)))*fslip(ifa)
      srd(:)=srtot(:)*ABS(SIN(fault_mech(ifa,3)))*fslip(ifa)
    ELSE IF (iread == 0.AND.dipole == 1) THEN
      srs(:)=srtot(:)*fslip(ifa)*0.5
      srd(:)=srtot(:)*fslip(ifa)*0.5
    END IF

    ! For slip scaling visualization
    IF (ifa.EQ.18503) THEN
      WRITE(70,*) srd(:)
      sltmp=0.0
      DO it=2,ntall
        tt=tt+dt
        sltmp=sltmp+dt*srd(it)
        WRITE(72,*) sltmp
      END DO
    ELSE IF (ifa.EQ.18534) THEN
      WRITE(71,*) srd(:)
      sltmp=0.0
      DO it=2,ntall
        tt=tt+dt
        sltmp=sltmp+dt*srd(it)
        WRITE(73,*) sltmp
      END DO
    END IF
!
! Current Strike and Dip Fault Angles
!
    str=fault_mech(ifa,1)
    dip=fault_mech(ifa,2)
!
! 6 Component Tensor (momentrate = sliprate*mu*area)
! Along-Strike Moment-Rate Tensor Components:
!
    rakes=0.0*pi/180.0
!
    ayy_s = -(SIN(dip)*COS(rakes)*SIN(2.*str)+ &
        SIN(2.*dip)*SIN(rakes)*SIN(str)*SIN(str))
!
    axy_s = SIN(dip)*COS(rakes)*COS(2.*str)+ &
        0.5*(SIN(2.*dip)*SIN(rakes)*SIN(2.*str))
!
    ayz_s = (COS(dip)*COS(rakes)*COS(str)+ &
        COS(2.*dip)*SIN(rakes)*SIN(str))
!
    axx_s = SIN(dip)*COS(rakes)*SIN(2.*str)- &
        SIN(2.*dip)*SIN(rakes)*COS(str)*COS(str)
!
    axz_s =  (COS(dip)*COS(rakes)*SIN(str)- &
        COS(2.*dip)*SIN(rakes)*COS(str))
!
    azz_s = SIN(2.*dip)*SIN(rakes)
!
    xx_s = axx_s*areamu(icont1)
    yy_s = ayy_s*areamu(icont1)
    zz_s = azz_s*areamu(icont1)
    xz_s = axz_s*areamu(icont1)
    yz_s = ayz_s*areamu(icont1)
    xy_s = axy_s*areamu(icont1)
!
! Along-Dip Moment-Rate Tensor Components:
!
    raked =90.0*pi/180.0
    !raked =270.0*pi/180.0
!
    ayy_d = -(SIN(dip)*COS(raked)*SIN(2.*str)+ &
        SIN(2.*dip)*SIN(raked)*SIN(str)*SIN(str))
!
    axy_d = SIN(dip)*COS(raked)*COS(2.*str)+ &
        0.5*(SIN(2.*dip)*SIN(raked)*SIN(2.*str))
!
    ayz_d = (COS(dip)*COS(raked)*COS(str)+ &
        COS(2.*dip)*SIN(raked)*SIN(str))
!
    axx_d = SIN(dip)*COS(raked)*SIN(2.*str)- &
        SIN(2.*dip)*SIN(raked)*COS(str)*COS(str)
!
    axz_d =  (COS(dip)*COS(raked)*SIN(str)- &
        COS(2.*dip)*SIN(raked)*COS(str))
!
    azz_d = SIN(2.*dip)*SIN(raked)
!
    xx_d = axx_d*areamu(icont1)
    yy_d = ayy_d*areamu(icont1)
    zz_d = azz_d*areamu(icont1)
    xz_d = axz_d*areamu(icont1)
    yz_d = ayz_d*areamu(icont1)
    xy_d = axy_d*areamu(icont1)
!
! Current Fault Point Coordinates 
!
    cx=(coord3d(ifa,1)-1)*dx
    cy=(coord3d(ifa,2)-1)*dx
    cz=(coord3d(ifa,3)-1)*dx
!
! Rupture Time of Current Fault Point (WARNING: Only in planar faults)
!
    dist=SQRT((cnx-cx)**2.0+(cny-cy)**2.0+(cnz-cz)**2.0)
    rtime=dist/rvel
    !WRITE(*,*) cz/1000.0,cx/1000.0,rtime
    WRITE(31,FMT=43) cz/1000.0,cx/1000.0,rtime
    IF (iread == 1) rtime=0.0
!
! Time Shift of Moment-Rate Functions Depending on Rupture Velocity
!
    IF (iout == 1) THEN
      IF (iread == 1) THEN
!       If read sliprate
        DO it=1,ntall
          WRITE(37,FMT=41) &
          xx_d*srd(it), &
          yy_d*srd(it), &
          zz_d*srd(it), &
          xz_d*srd(it), &
          yz_d*srd(it), &
          xy_d*srd(it)
!         xx_s*srs(it) + xx_d*srd(it), &
!         yy_s*srs(it) + yy_d*srd(it), &
!         zz_s*srs(it) + zz_d*srd(it), &
!         xz_s*srs(it) + xz_d*srd(it), &
!         yz_s*srs(it) + yz_d*srd(it), &
!         xy_s*srs(it) + xy_d*srd(it)
        END DO
      ELSE
!       If no read sliprate
        tt=0.0
        icont=1
        IF (dipole == 0) THEN
          DO it=1,ntall
            IF (tt >= rtime) THEN
              WRITE(37,FMT=41) &
              xx_s*srs(icont) + xx_d*srd(icont), &
              yy_s*srs(icont) + yy_d*srd(icont), &
              zz_s*srs(icont) + zz_d*srd(icont), &
              xz_s*srs(icont) + xz_d*srd(icont), &
              yz_s*srs(icont) + yz_d*srd(icont), &
              xy_s*srs(icont) + xy_d*srd(icont)
              icont=icont+1
            ELSE
              WRITE(37,FMT=41) 0.0,0.0,0.0,0.0,0.0,0.0
            END IF
            tt=tt+dt
          END DO
        ELSE
          DO it=1,ntall
            IF (tt >= rtime) THEN
              WRITE(37,FMT=41) &
              0.0, &
              0.0, &
              areamu*srs(icont) + areamuctte*srd(icont), &
              0.0, &
              0.0, &
              0.0
              icont=icont+1
            ELSE
              WRITE(37,FMT=41) 0.0,0.0,0.0,0.0,0.0,0.0
            END IF
            tt=tt+dt
          END DO
        END IF
      END IF
    END IF
  END DO
!
  PRINT *
  WRITE(*,*) 'END PROGRAM'
  PRINT *
  WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT *
!
  40 FORMAT(3I5)
  41 FORMAT(6E15.7)
  42 FORMAT(2F12.5)
  43 FORMAT(3F12.5)
!
  CONTAINS
!
! ++++++++++++++++++++++++++++++
!
  FUNCTION trian(tr,tt,t0,dt,pi)
!
  IMPLICIT NONE
!
  REAL :: tr,tt,t0,dt
  REAL :: a0,pi,trian
!
  t0=0.0
  a0=2.0/tr
  IF ((tt-t0) <= tr/2.0) THEN
    trian=(2.0*a0/tr)*(tt-t0)
  ELSE IF ((tt-t0) > tr/2.0 .AND. (tt-t0) <= tr) THEN
    trian=-(2.0*a0/tr)*(tt-t0)+2.0*a0
  ELSE IF ((tt-t0) > tr) THEN
    trian=0.0
  END IF
!
  END FUNCTION trian
!
! ++++++++++++++++++++++++++++++
!
  FUNCTION gauss(tr,tt,t0,dt,pi)
!
  IMPLICIT NONE
!
  REAL :: tr,tr4,tt,t0,dt
  REAL :: pi,gauss
!
  tr4=tr/4.0
  t0=2.0*tr
! t0=tr
  gauss=1.0/tr4/SQRT(pi)*EXP(-(tt-t0)**2.0/tr4**2.0)
!
  END FUNCTION gauss
!
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
END PROGRAM
