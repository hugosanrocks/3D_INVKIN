!

SUBROUTINE subcell(xf0,yf0,zf0,ang,xlf,h,nsp,coord2d,ifpoint, &
           ipoint,iseg,pentx,ordx,xo,zo,xf,zf,m,dist,ind,indxo, &
           indzo,fseg,coord3d,vnormal,ifp3d,iaxis,wide,nlength, &
           nwide,nori,indyo,fault_mech,rake)
!
  IMPLICIT NONE
!
! Variable declaration
!
  INTEGER :: i,j,k
  INTEGER :: indxo,indzo,indyo
  INTEGER :: ipoint,ifpoint,inz,ifp3d
  INTEGER :: indx,indz
  INTEGER :: nsp
  INTEGER :: iseg,fseg,node
  INTEGER :: ind(3,2)
  INTEGER :: coord2d(ifpoint,2)
  INTEGER :: coord3d(ifpoint,3)
  INTEGER :: isignx,isignz,ico1,ico2
  INTEGER :: iaxis,nlength,nwide
  INTEGER :: icont2,nori,idist
  REAL :: ang,h
  REAL :: xf0,yf0,zf0,xlf
  REAL :: pi,radi,radia
  REAL :: pentx,pentz,MIN
  REAL :: ordx,ordz
  REAL :: pent,orde
  REAL :: xo,zo,xf,zf
  REAL :: xot,zot
  REAL :: delf,fin,nfin
  REAL :: m(4),dist(3)
  REAL :: vnormal(ifpoint,3)
  REAL :: fault_mech(ifpoint,3)
  REAL :: wide,rake
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Set constants
!
  pi=3.141592654
  radi=ang*pi/180.0
  radia=radi
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! The initial node
!
  IF (iseg == 1) THEN
!
    xot=INT(xf0/h)*h
    zot=INT(zf0/h)*h
!
    m(1)=SQRT((xf0-xot)**2.0+(zf0-zot)**2.0)
    m(2)=SQRT((xf0-xot-h)**2.0+(zf0-zot)**2.0)
    m(3)=SQRT((xf0-xot-h)**2.0+(zf0-zot-h)**2.0)
    m(4)=SQRT((xf0-xot)**2.0+(zf0-zot-h)**2.0)
!
    MIN=10E4
    DO i=1,4
      IF (m(i) < MIN) MIN=m(i)
    END DO
!
    IF (m(1) == MIN) THEN
      xo=xot
      zo=zot
    ELSE IF (m(2) == MIN) THEN
      xo=xot+h
      zo=zot
    ELSE IF (m(3) == MIN) THEN
      xo=xot+h
      zo=zot+h
    ELSE IF (m(4) == MIN) THEN
      xo=xot
      zo=zot+h
    END IF
!
    indxo=NINT(xo/h)+1
    indzo=NINT(zo/h)+1
!
    ipoint=ipoint+1
!   coord2d(ipoint,1)=indxo+nsp
!   coord2d(ipoint,2)=indzo+nsp
    coord2d(ipoint,1)=indxo
    coord2d(ipoint,2)=indzo
!
    xo=indxo*h
    zo=indzo*h
    GO TO 100
!
  END IF
!
! Set line parameters
!
  xo=xf
  zo=zf
!
  100     CONTINUE
!
  IF (iaxis == 3) THEN
    xf=xo+xlf*COS(radi)
    zf=zo+xlf*SIN(radi)
  ELSE
    xf=xo+xlf*COS(radi)
    zf=zo-xlf*SIN(radi)
  END IF
!
  pentx=(zf-zo)/(xf-xo)
  pentz=(xf-xo)/(zf-zo)
  ordx=zo-pentx*xo
  ordz=xo-pentz*zo
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Searching the following nodes
!
  IF (ang >= 0.0.AND.ang < 90.0) THEN
    IF (iaxis == 3) THEN
      isignx=1
      isignz=1
    ELSE
      isignx=1
      isignz=-1
    END IF
  END IF
  IF (ang >= 90.0.AND.ang < 180.0) THEN
    IF (iaxis == 3) THEN
      isignx=-1
      isignz=1
    ELSE
      isignx=-1
      isignz=-1
    END IF
  END IF
  IF (ang >= 180.0.AND.ang < 270.0) THEN
    IF (iaxis == 3) THEN
      isignx=-1
      isignz=-1
    ELSE
      isignx=-1
      isignz=1
    END IF
  END IF
  IF (ang >= 270.AND.ang <= 360.0) THEN
    IF (iaxis == 3) THEN
      isignx=1
      isignz=-1
    ELSE
      isignx=1
      isignz=1
    END IF
  END IF
!
  IF ((ang > 45.0.AND.ang <= 135.0).OR. &
      (ang >= 225.0.AND.ang <= 315.0)) THEN
    pent=pentz
    orde=ordz
    fin=zf
    ico1=1
    ico2=2
    radi=ABS(radi-pi/2.0)
!
  ELSE IF ((ang >= 0.0.AND.ang <= 45.0).OR.                             &
        (ang > 135.0.AND.ang < 225.0).OR.                               &
        (ang > 315.0.AND.ang <= 360.0)) THEN
    pent=pentx
    orde=ordx
    fin=xf
    ico1=2
    ico2=1
  END IF
!
  DO j=1,ifpoint
!
    indx=indxo
    indz=indzo
!
    ind(1,1)=indx
    ind(1,2)=indz+isignz
    dist(1)=ABS((h*ind(1,ico1)-(h*ind(1,ico2)* &
            pent+orde))*ABS(COS(radi)))
!
    ind(2,1)=indx+isignx
    ind(2,2)=indz+isignz
    dist(2)=ABS((h*ind(2,ico1)-(h*ind(2,ico2)* &
            pent+orde))*ABS(COS(radi)))
!
    ind(3,1)=indx+isignx
    ind(3,2)=indz
    dist(3)=ABS((h*ind(3,ico1)-(h*ind(3,ico2)* &
            pent+orde))*ABS(COS(radi)))
!
    MIN=10E4
    DO i=1,3
      IF (dist(i) < MIN) THEN
        MIN=dist(i)
        idist=i
      END IF
    END DO
!
    DO i=1,3
      IF (idist == i) THEN
        indxo=ind(i,1)
        indzo=ind(i,2)
        ipoint=ipoint+1
        coord2d(ipoint,1)=indxo
        coord2d(ipoint,2)=indzo
!       coord2d(ipoint,1)=indxo+nsp
!       coord2d(ipoint,2)=indzo+nsp
        node=i
        GO TO 300
      END IF
    END DO
    300        CONTINUE
    IF (xlf < h) THEN
      ipoint=ipoint-1
    END iF
!
    delf=ABS(fin-ind(node,ico2)*h)
    nfin=ind(node,ico2)*h
!
    IF (iseg == fseg) THEN
      IF (delf < h) GO TO 200
    ELSE
      IF (ico1 == 1) THEN
        IF (delf <= h.AND.nfin > fin) GO TO 200
      ELSE
        IF (delf <= h.AND.nfin <= fin) GO TO 200
      END IF
    END IF
!
  END DO
  200     CONTINUE
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Create the third dimension
!
! indyo=NINT(yf0/h)+nsp
  indyo=NINT(yf0/h)+1
  nwide=NINT(wide/h)+1
!
  icont2=0
  DO i=nori,ipoint
    DO j=1,nwide
      ifp3d=ifp3d+1
      IF (iaxis == 1) THEN
        coord3d(ifp3d,1)=indyo+icont2
        coord3d(ifp3d,2)=coord2d(i,1)
        coord3d(ifp3d,3)=coord2d(i,2)
!
        vnormal(ifp3d,1)=0.0
        vnormal(ifp3d,2)=SIN(radia)
        vnormal(ifp3d,3)=COS(radia)
!
!       Fault Strike
        IF (ang > 0.0.AND.ang < 90.0.OR. &
            ang > 180.AND.ang < 270.0) THEN
          fault_mech(ifp3d,1)=270.0
        ELSE
          fault_mech(ifp3d,1)=90.0
        END IF
!
!       Fault Dip
        IF (ang >= 0.0.AND.ang < 90.0) fault_mech(ifp3d,2)=ang
        IF (ang >= 90.0.AND.ang < 180.0) fault_mech(ifp3d,2)=180.0-ang
        IF (ang >= 180.0.AND.ang < 270.0) fault_mech(ifp3d,2)=ang-180.0
        IF (ang >= 270.0.AND.ang <= 360.0) fault_mech(ifp3d,2)=360.0-ang
!
      END IF
      IF (iaxis == 2) THEN
        coord3d(ifp3d,1)=coord2d(i,1)
        coord3d(ifp3d,2)=indyo+icont2
        coord3d(ifp3d,3)=coord2d(i,2)
!
        vnormal(ifp3d,1)=SIN(radia)
        vnormal(ifp3d,2)=0.0
        vnormal(ifp3d,3)=COS(radia)
!
!       Fault Strike
        IF (ang > 0.0.AND.ang < 90.0.OR. &
            ang >= 180.AND.ang < 270.0) THEN
          fault_mech(ifp3d,1)=180.0
        ELSE
          fault_mech(ifp3d,1)=0.0
        END IF
!
!       Fault Dip
        IF (ang >= 0.0.AND.ang < 90.0) fault_mech(ifp3d,2)=ang
        IF (ang >= 90.0.AND.ang < 180.0) fault_mech(ifp3d,2)=180.0-ang
        IF (ang >= 180.0.AND.ang < 270.0) fault_mech(ifp3d,2)=ang-180.0
        IF (ang >= 270.0.AND.ang <= 360.0) fault_mech(ifp3d,2)=360.0-ang
!
      END IF
      IF (iaxis == 3) THEN
        coord3d(ifp3d,1)=coord2d(i,1)
        coord3d(ifp3d,2)=coord2d(i,2)
        coord3d(ifp3d,3)=indyo+icont2
!
        vnormal(ifp3d,3)=0.0
        vnormal(ifp3d,1)=SIN(radia)
        vnormal(ifp3d,2)=-COS(radia)
!
!       Fault Strike
        IF (ang >= 0.0.AND.ang < 90.0) THEN
          fault_mech(ifp3d,1)=90.0-ang
        ELSE
          fault_mech(ifp3d,1)=360.0+90.0-ang
        END IF
!
!       Fault Dip
        fault_mech(ifp3d,2)=90.0
      END IF
!     Fault Dip
      fault_mech(ifp3d,3)=rake
      icont2=icont2+1
    END DO
    icont2=0
  END DO
!
  nlength=ipoint
  nori=ipoint+1
!
  WRITE(*,20)iseg,xo,zo
  WRITE(*,30)iseg,xf,zf
!
  20 FORMAT(' Origin of segment',i3,' at Xo = ',f8.1,' and Zo = ',f8.1)
  30 FORMAT('    End of segment',i3,' at Xf = ',f8.1,' and Zf = ',f8.1)
  40 FORMAT(3I4)
!
  RETURN
END SUBROUTINE subcell
