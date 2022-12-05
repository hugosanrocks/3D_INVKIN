************************************************************************
* convertion geographique (lon,lat) en Lambert(x,y) pour les coordonnees
* des stations
************************************************************************

      PROGRAM pre_code

      IMPLICIT NONE

      INTEGER k,nb
      INTEGER nstk,ndip,idip,istk,is,ns
      REAL xlat0,xlon0,Y,X,Z,xlat,xlon
      REAL t(100000),xs(100000),ys(100000),zs(100000)
      REAL pi
      REAL astk,adip,estk,edip
      REAL phi,delta,rake,phi_axi
      REAL slip
      REAL dstk,ddip
      REAL dxs,dxd,dys,dyd,dzd
      REAL zhypo
      REAL vr
      CHARACTER*4 name

      OPEN(10,file='stations.in')
      OPEN(11,file='fault.in')
      OPEN(12,file='stations.out')
      OPEN(13,file='stations.name')
      OPEN(14,form='formatted',file='fault_3d.out')
      OPEN(15,form='formatted',file='stations_2d.out')
c     OPEN(15,form='formatted',file='axi.hist')
      OPEN(16,form='formatted',file='fault_2d.out')

      pi=3.14159265359

c +++++++++++++++++++++++++++++++
c Read and process stations data
c +++++++++++++++++++++++++++++++

      READ(10,*) xlat0,xlon0
      WRITE(6,*) 'latitude longitude de l''epicentre'
      WRITE(6,*) xlat0,xlon0


      READ(10,*) nb

      !WRITE(13,*) 'EPIC',xlat0,xlon0
      WRITE(6,*) 'nombre de stations'
      WRITE(6,*) nb
      WRITE(6,*) 'latitude, longitude, altitude(m), nom(4 lettres)'
      WRITE(6,*) 'x,y,lon,lat'

      CALL InitGlobalLambert(xlat0,xlon0)
      DO k=1,nb
         READ(10,*) xlat,xlon,Z,name
         !CALL GeoToGlobalLambert(xlat,xlon,Y,X)
         CALL GeoToGlobalLambert(xlat,xlon,X,Y)
         WRITE(6,*) name,X,Y,xlon,xlat
         WRITE(12,*) 1000*X,1000*Y,1000*Z 
         WRITE(15,*) 1000*X,1000*Y
         WRITE(13,'(a4)') name 
      END DO
      CLOSE(12)
      CLOSE(13)

c +++++++++++++++++++++++++++++++
c Read and process fault data
c +++++++++++++++++++++++++++++++

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) 'MAKe Sure The SOURCE file has been deleted first'
      WRITE(*,*)
      WRITE(*,*)
       
      READ(11,*) astk,adip
      WRITE(6,*) 'longueur(stk) largeur(dip) en km        ',astk,adip
      WRITE(*,*) astk,adip
      READ(11,*) estk, edip
      WRITE(6,*) 'position de l''epic. sur la faille en km ',estk,edip
      WRITE(*,*)estk, edip
      READ(11,*) zhypo
      WRITE(6,*)"profondeur de l'hypocentre en km:         ", zhypo
      WRITE(*,*) zhypo
      READ(11,*) phi,delta,rake
      WRITE(6,*)"azim faille en degre par rapport au nord ",phi
      WRITE(*,*)phi
      WRITE(6,*)"pendage (dip) en degre:                  ",delta
      WRITE(*,*) delta
      WRITE(6,*)"rake en degre:                           ",rake
      WRITE(*,*) rake
      READ(11,*) nstk,ndip
      WRITE(6,*)"nb de sousfailles selon strk et dip:     ",nstk,ndip
      WRITE(*,*) nstk,ndip
c      READ(11,*) slip
c      WRITE(6,*)"amplitude du glissement (m) ?            ", slip
c      WRITE(*,*) slip
c      READ(11,*) vr
c      WRITE(6,*)"vitesse de rupture en m/s:               ", vr
c      WRITE(*,*)vr
c      WRITE(*,*)

      CLOSE(11)
      phi_axi = 270-phi
      estk = astk-estk 
      astk = astk*1e3
      edip = edip*1e3
      estk = estk*1e3
      zhypo = zhypo*1e3
      adip = adip*1e3
************************************************************************
* nombre de points source (indexes de 1 a ns)
************************************************************************

      ns=nstk*ndip

***********************************************************************
* Modified by John Diaz to adjust to DWN orientation
***********************************************************************

      phi = 270-phi

!     IF (rake .lt. 90)THEN
!       rake = 120+rake
!       ELSEIF(rake .gt. 90 .and. rake .le. 180)THEN
!             rake = 180-rake
!       ELSEIF(rake .gt. 180)THEN
!           rake = 540-rake
!     END IF 

! Fixed by Victor
      IF (rake .le. 180)THEN
        rake = 180-rake
      ELSE
        rake = 540-rake
      END IF 

!      print*, ' phi ',phi,' rake ', rake  
************************************************************************
* angles en radian
************************************************************************

      phi=phi*pi/180.
      delta=delta*pi/180.

************************************************************************
* taille de chaque sous faille selon le strike et le dip
************************************************************************

      dstk=astk/nstk
      ddip=adip/ndip
      WRITE(6,*)"dstk, dsty = ",dstk,ddip

************************************************************************
* decalage selon x,y ou z provoque entre sous-failles contigues
************************************************************************

      dxs=dstk*cos(phi)
      dxd=-ddip*cos(delta)*sin(phi)
      dys=dstk*sin(phi)
      dyd=ddip*cos(delta)*cos(phi)
      dzd=ddip*sin(delta)
c      WRITE(6,*)"5"

************************************************************************
* coordonnees de la premiere sous-faille dans le repere x,y et z
************************************************************************

      xs(1)=(dxs+dxd)/2-estk*cos(phi)+edip*cos(delta)*sin(phi)
      ys(1)=(dys+dyd)/2-estk*sin(phi)-edip*cos(delta)*cos(phi)
      zs(1)=dzd/2-edip*sin(delta)+zhypo
      WRITE(6,*)"first = (",xs(1), ys(1), zs(1)," )"

************************************************************************
* coordonnees des autres sous-failles et ecriture dans les fichiers
************************************************************************

      is=1
      do idip=1,ndip
        do istk=1,nstk
c        WRITE(6,*)"1"
          xs(is)=xs(1)+(istk-1)*dxs+(idip-1)*dxd
          ys(is)=ys(1)+(istk-1)*dys+(idip-1)*dyd
          zs(is)=zs(1)+(idip-1)*dzd
c         t(is)=(((xs(is)**2)+(ys(is)**2)+((zs(is)-zhypo)**2))**0.5)/vr
          WRITE(14,*) is,xs(is),ys(is),zs(is)
          WRITE(16,*) xs(is),ys(is)
          t(is)=t(is)
c         WRITE(15,'(I5,7G12.3)')is,slip,phi_axi,delta*180/pi,
c    &                           rake,ddip,dstk,t(is)
          is=is+1
       enddo
      enddo
      WRITE(6,*)"last = (",xs(is-1), ys(is-1), zs(is-1)," )"


      CLOSE(14)
      CLOSE(15)
      CLOSE(16)

************************************************************************

      END


      SUBROUTINE InitGlobalLambert(lat0,lon0)

c       PROJECTION LAMBERT CALCULEE EN SPHERIQUE
c       AVEC LAT0 ET LON0 COMME CENTRE DE PROJECTION }

      implicit none

      REAL lat0,lon0

      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,initialized
      logical initialized
      REAL L0,latitude0,longitude0,n,pi,r0

        if (lon0.lt.0) lon0=lon0+360
      if (abs(lon0).ge.360.or.abs(lat0).ge.90) then
        print *,'InitGlobalLambert: parameter values out of range'
        stop
      endif
      pi=atan(1.)*4
      latitude0=lat0*pi/180
      longitude0=lon0*pi/180
      if (abs(latitude0).lt.1e-6) latitude0=1e-6
      if (abs(latitude0-pi/2).lt.1e-6) latitude0=pi/2-1e-6
      if (abs(latitude0+pi/2).lt.1e-6) latitude0=-pi/2+1e-6
      n=sin(latitude0)
      r0=6371/tan(latitude0)
      L0=log(tan(pi/4+latitude0/2))
      initialized=.true.
      end
!
!  SUBROUTINE GeoToGlobalLambert(la,lo,X,Y)
!
      SUBROUTINE GeoToGlobalLambert(la,lo,X,Y)
      
      implicit none

      REAL X,Y,la,lo

      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,
     &         initialized
      logical initialized
      REAL L0,latitude0,longitude0,n,pi,r0
      
      REAL L,gamma,r,lat,lon

      lat=la
      lon=lo
      if (.not.initialized) then
        print *,'GeoToGlobalLambert: InitGlobalLambert missing'
        stop
      endif
      if (lon.lt.0 ) lon=lon+360
      if (abs(lon).ge.360.or.abs(lat).ge.90) then
        print *,'GeoToGlobalLambert: parameter values out of range'
        stop
      endif
      lon=lon*pi/180
      lat=lat*pi/180
      L=log(tan(pi/4+lat/2))
      gamma=n*(lon-longitude0)
      r=r0*exp(n*(L0-L))
      X=r*sin(gamma)
      Y=R0-r*cos(gamma)
      end
c      { -------------------- }
      SUBROUTINE GlobalLambertToGeo(X,Y,lat,lon)

      implicit none

      REAL X,Y,lat,lon
      
      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,
     &                        initialized
      logical initialized
      REAL L0,latitude0,longitude0,n,pi,r0

      REAL L,gamma,r

      if (.not.initialized) then
        print *,'GlobalLambertToGeo: InitGlobalLambert missing'
        stop
      endif
      if (abs(r0-Y).lt.1e-6) then
        gamma=pi/2
        if (X.lt.0) gamma=-gamma
      else
        gamma=atan(X/(r0-Y))
        if (abs(gamma).lt.1e-6) then 
          r=r0-Y 
        else 
          r=X/sin(gamma)
        endif
      endif
      lat=atan(exp(L0-log(r/r0)/n))*360/pi-90
      lon=(longitude0+gamma/n)*180/pi
      if (lon.lt.-180) then 
        lon=lon+360 
      else if (lon.gt.180) then
        lon=lon-360
      endif
      end
c      { -------------------- }
