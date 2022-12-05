!

SUBROUTINE geometry(ifpoint,ifp3d,inp3d,h, &
                    coord3d,vnormal,nucle, &
                    fault_mech,ixs,iys,izs, &
                    nlength,nwide,coord2d)
!
  IMPLICIT NONE
!
  INTEGER :: i,j,k
  INTEGER :: nsp
  INTEGER :: ixs,iys,izs
  INTEGER :: ifpoint
  INTEGER :: ipoint
  INTEGER :: ifp3d,inp3d
  INTEGER :: nlength,nwide,nori
  INTEGER :: iaxis
  INTEGER :: coord2d(ifpoint,2)
  INTEGER :: coord3d(ifpoint,3)
  INTEGER :: nucle(ifpoint,4)
! INTEGER :: inde(3,2)
  INTEGER :: inde(7,2)
  INTEGER :: ijoin,ibreak
  INTEGER :: itmp,iseg
  INTEGER :: indxo,indzo,indyo
  INTEGER :: ncn1,ncn2
  INTEGER :: npoint
  REAL :: vnormal(ifpoint,3)
  REAL :: fault_mech(ifpoint,3)
! REAL :: m(4),dist(3)
  REAL :: m(4),dist(7)
  REAL :: pentx,ordx,xo,zo,xf,zf
  REAL :: ang,h
  REAL :: cn1,cn2
  REAL :: flong,xf0,zf0,yf0
  REAL :: wide,rake
  REAL :: radius,radi2
  REAL :: hdist,vdist
  REAL :: xdist,ydist,zdist
  REAL :: hcell
!
  PARAMETER (nsp=15)      	! Not active (PML)
!
! Nucleation zone
!
  DO i=1,6
    READ(33,*)
  END DO
  READ(33,*) cn1,cn2
  cn1=cn1*1e3
  cn2=cn2*1e3
  REWIND(33)
!
  hcell=h 
  ipoint=0
  nori=1
  ifp3d=0
  inp3d=0
!
  coord2d(:, :)=0
  coord3d(:, :)=0
  vnormal(:, :)=0.0
  fault_mech(:, :)=0.0
  dist(:)=0.0
  inde(:, :)=0.0
  inde(:, :)=0.0
  m(:)=0.0
  zf=0.0
  xf=0.0
!
  DO i=1,4
    READ(33,*)
  END DO
  READ(33,*)ibreak
  READ(33,*)iaxis
  DO i=1,6
    READ(33,*)
  END DO
!
  DO j=1,ibreak
    READ(33,*)ijoin
    IF (iaxis == 1) READ(33,*) yf0,xf0,zf0
    IF (iaxis == 2) READ(33,*) xf0,yf0,zf0
    IF (iaxis == 3) READ(33,*) xf0,zf0,yf0
    READ(33,*)wide
    xf0=xf0*1e3
    yf0=yf0*1e3
    zf0=zf0*1e3
    wide=wide*1e3
    itmp=j
    iseg=1
    DO i=1,ijoin
      IF (j == itmp) THEN
        WRITE(*,*)
        WRITE(*,40) j
        PRINT *
      END IF
      READ(33,*) ang,flong,rake
      flong=flong*1e3
!
      CALL subcell(xf0,yf0,zf0,ang,flong,h,nsp,coord2d, &
                   ifpoint,ipoint,iseg,pentx,ordx,xo,zo, &
                   xf,zf,m,dist,inde,indxo,indzo,ijoin, &
                   coord3d,vnormal,ifp3d,iaxis,wide, &
                   nlength,nwide,nori,indyo,fault_mech,rake)
      iseg=iseg+1
      itmp=0
    END DO
    IF (j.LT.ibreak) THEN
      DO i=1,5
        READ(33,*)
      END DO
    END IF
  END DO
!
! Create square nucleation zone
!
  ncn1=NINT(cn1/hcell)+1
  ncn2=NINT(cn2/hcell)+1
  npoint=nwide*(ncn1-1)+ncn2
!
  PRINT *
  WRITE(*,*) 'Nucleation Point Coordinates:'
  WRITE(*,FMT=44) (ncn1-1)*h
  WRITE(*,FMT=45) (ncn2-1)*h
!
  ixs=coord3d(npoint,1)
  iys=coord3d(npoint,2)
  izs=coord3d(npoint,3)
!
  radius=2.0*h     	! Nucleation Radius
  radi2=(radius/2.0)-(hcell/2.0)
  DO i=1,ifp3d
    xdist=ABS(coord3d(i,1)-ixs)*h
    ydist=ABS(coord3d(i,2)-iys)*h
    zdist=ABS(coord3d(i,3)-izs)*h
    IF (xdist <= radi2 .AND. ydist <= radi2 .AND. zdist <= radi2) THEN
      inp3d=inp3d+1
      nucle(inp3d,1)=coord3d(i,1)
      nucle(inp3d,2)=coord3d(i,2)
      nucle(inp3d,3)=coord3d(i,3)
      nucle(inp3d,4)=i
    END IF
  END DO
!
  WRITE(*,*)
  WRITE(*,*)'Amount of fault points =',ifp3d
  WRITE(*,*)'Amount of nucleation points =',inp3d
  WRITE(*,*)
  WRITE(*,FMT=47) nlength
  WRITE(*,FMT=46) nwide
  PRINT *
!
! Write out fault coordinates
!
  DO i=1,nlength
    WRITE(34,FMT=41) coord2d(i,1),coord2d(i,2)
  END DO
!
  DO i=1,ifp3d
    WRITE(35,FMT=42) coord3d(i,1),coord3d(i,2),coord3d(i,3)
    WRITE(36,FMT=43) fault_mech(i,1),fault_mech(i,2),fault_mech(i,3)
  END DO
!
  40 FORMAT (' Coordinates of fault number',i2,':')
  41 FORMAT (2I5)
  42 FORMAT (3I5)
  43 FORMAT (3F12.5)
  44 FORMAT ('   Along the non-planar direction: ',F10.2,' m')
  45 FORMAT ('   Along the translation invariant direction: ',F10.2,' m')
  46 FORMAT (' Fault points in translation invariant direction =',I6)
  47 FORMAT ('            Fault points in non-planar direction =',I6)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  RETURN
END SUBROUTINE geometry
