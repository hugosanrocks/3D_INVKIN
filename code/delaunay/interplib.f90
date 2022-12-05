

      SUBROUTINE FIND1D(Xwant,Value,X,Array,NX)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: NX
      REAL :: Value , Xwant
      REAL , DIMENSION(NX) :: Array , X
      INTENT (IN) Array , NX
      INTENT (OUT) Value
!
! Local variables
!
      REAL :: c, m
      INTEGER :: jx , jxp1
!
!-----------------------------------------------------------------------
!
      CALL POS(X,NX,Xwant,jx)
      IF ( jx==0 ) THEN
         jx = 1
      ELSEIF ( jx>=NX ) THEN
         jx = NX - 1
      ENDIF
      jxp1 = jx + 1
!
      m = (Array(jx+1)-Array(jx))/(X(jx+1)-X(jx))
      c = Array(jx) - (m*X(jx))
      !print *, m, c, array(jx+1), array(jx), x(jx+1), x(jx), xwant
      Value = Xwant*m + c
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIND1D
      SUBROUTINE FIND2D(Xwant,Ywant,Value,X,Y,Array,NX,NY)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: NX , NY
      REAL :: Value , Xwant , Ywant
      REAL , DIMENSION(NX,NY) :: Array
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      INTENT (IN) NX , NY
!
! Local variables
!
      INTEGER :: jx , jxp1 , jy , jyp1
!
!-----------------------------------------------------------------------
!
!     written by:   David Lary
!
!     started:      7/1/1993
!
!     last updated: 22/1/2004
!
!----------------------------------------------------------------------
!
!     2D bilinear interpolation.
!
!-----------------------------------------------------------------------
!
      CALL FINDINDEX(Xwant,Ywant,X,Y,NX,NY,jx,jxp1,jy,jyp1)
      CALL FINDPOINT(Xwant,Ywant,Value,X,Y,Array,NX,NY,jx,jxp1,jy,jyp1)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIND2D
      SUBROUTINE FIND3D(Xwant,Ywant,Zwant,Value,X,Y,Z,Array,NX,NY,NZ)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: NX , NY , NZ
      REAL :: Value , Xwant , Ywant , Zwant
      REAL , DIMENSION(NX,NY,NZ) :: Array
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) NX , NY , NZ
!
! Local variables
!
      INTEGER :: jx , jxp1 , jy , jyp1 , jz , jzp1
!
!----------------------------------------------------------------------
!
      CALL FINDINDEX3D(Xwant,Ywant,Zwant,X,Y,Z,NX,NY,NZ,jx,jxp1,jy,jyp1,&
                     & jz,jzp1)
!
      CALL FINDPOINT3D(Xwant,Ywant,Zwant,Value,X,Y,Z,Array,NX,NY,NZ,jx, &
                     & jxp1,jy,jyp1,jz,jzp1)
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FIND3D
      SUBROUTINE FIND4D(Xwant,Ywant,Zwant,Qwant,Value,X,Y,Z,Q,Array,NX, &
                      & NY,NZ,NQ)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: NQ , NX , NY , NZ
      REAL :: Qwant , Value , Xwant , Ywant , Zwant
      REAL , DIMENSION(NX,NY,NZ,NQ) :: Array
      REAL , DIMENSION(NQ) :: Q
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) NQ , NX , NY , NZ
!
! Local variables
!
      INTEGER :: jq , jqp1 , jx , jxp1 , jy , jyp1 , jz , jzp1
!
!----------------------------------------------------------------------
!
      CALL FINDINDEX4D(Xwant,Ywant,Zwant,Qwant,X,Y,Z,Q,NX,NY,NZ,NQ,jx,  &
                     & jxp1,jy,jyp1,jz,jzp1,jq,jqp1)
!
      CALL FINDPOINT4D(Xwant,Ywant,Zwant,Qwant,Value,X,Y,Z,Q,Array,NX,  &
                     & NY,NZ,NQ,jx,jxp1,jy,jyp1,jz,jzp1,jq,jqp1)
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FIND4D
      SUBROUTINE FINDINDEX(Xwant,Ywant,X,Y,NX,NY,Jx,Jxp1,Jy,Jyp1)
      IMPLICIT NONE
!
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1
      INTEGER :: NX , NY
      REAL :: Xwant , Ywant
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      INTENT (IN) NX , NY
      INTENT (OUT) Jxp1 , Jyp1
      INTENT (INOUT) Jx , Jy
!
!-----------------------------------------------------------------------
!
!     written by:   David Lary
!
!     started:      7/1/1993
!
!     last updated: 22/1/2004
!
!----------------------------------------------------------------------
!
!     Find grid references.
!
!-----------------------------------------------------------------------
!
!     i) x
      CALL POS(X,NX,Xwant,Jx)
      IF ( Jx==0 ) THEN
         Jx = 1
      ELSEIF ( Jx>=NX ) THEN
         Jx = NX - 1
      ENDIF
      Jxp1 = Jx + 1
!
!     ii) y
      CALL POS(Y,NY,Ywant,Jy)
      IF ( Jy==0 ) THEN
         Jy = 1
      ELSEIF ( Jy>=NY ) THEN
         Jy = NY - 1
      ENDIF
      Jyp1 = Jy + 1
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FINDINDEX
      SUBROUTINE FINDINDEX2D(Xwant,Ywant,X,Y,NX,NY,Jx,Jxp1,Jy,Jyp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1
      INTEGER :: NX , NY
      REAL :: Xwant , Ywant
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      INTENT (IN) NX , NY
      INTENT (OUT) Jxp1 , Jyp1
      INTENT (INOUT) Jx , Jy
!
!-----------------------------------------------------------------------
!
!     Find grid references.
!     i) x
      CALL POS(X,NX,Xwant,Jx)
      IF ( Jx==0 ) THEN
         Jx = 1
      ELSEIF ( Jx>=NX ) THEN
         Jx = NX - 1
      ENDIF
      Jxp1 = Jx + 1
!
!     ii) y
      CALL POS(Y,NY,Ywant,Jy)
      IF ( Jy==0 ) THEN
         Jy = 1
      ELSEIF ( Jy>=NY ) THEN
         Jy = NY - 1
      ENDIF
      Jyp1 = Jy + 1
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FINDINDEX2D
      SUBROUTINE FINDINDEX3D(Xwant,Ywant,Zwant,X,Y,Z,NX,NY,NZ,Jx,Jxp1,  &
                           & Jy,Jyp1,Jz,Jzp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1 , Jz , Jzp1
      INTEGER :: NX , NY , NZ
      REAL :: Xwant , Ywant , Zwant
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) NX , NY , NZ
      INTENT (OUT) Jxp1 , Jyp1 , Jzp1
      INTENT (INOUT) Jx , Jy , Jz
!
!----------------------------------------------------------------------
!
      CALL POS(X,NX,Xwant,Jx)
      IF ( Jx==0 ) THEN
         Jx = 1
      ELSEIF ( Jx>=NX ) THEN
         Jx = NX - 1
      ENDIF
      Jxp1 = Jx + 1
!
      CALL POS(Y,NY,Ywant,Jy)
      IF ( Jy==0 ) THEN
         Jy = 1
      ELSEIF ( Jy>=NY ) THEN
         Jy = NY - 1
      ENDIF
      Jyp1 = Jy + 1
!
      CALL POS(Z,NZ,Zwant,Jz)
      IF ( Jz==0 ) THEN
         Jz = 1
      ELSEIF ( Jz>=NZ ) THEN
         Jz = NZ - 1
      ENDIF
      Jzp1 = Jz + 1
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FINDINDEX3D
      SUBROUTINE FINDINDEX4D(Xwant,Ywant,Zwant,Qwant,X,Y,Z,Q,NX,NY,NZ,  &
                           & NQ,Jx,Jxp1,Jy,Jyp1,Jz,Jzp1,Jq,Jqp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jq , Jqp1 , Jx , Jxp1 , Jy , Jyp1 , Jz , Jzp1
      INTEGER :: NQ , NX , NY , NZ
      REAL :: Qwant , Xwant , Ywant , Zwant
      REAL , DIMENSION(NQ) :: Q
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) NQ , NX , NY , NZ
      INTENT (OUT) Jqp1 , Jxp1 , Jyp1 , Jzp1
      INTENT (INOUT) Jq , Jx , Jy , Jz
!
!----------------------------------------------------------------------
!
      CALL POS(X,NX,Xwant,Jx)
      IF ( Jx==0 ) THEN
         Jx = 1
      ELSEIF ( Jx>=NX ) THEN
         Jx = NX - 1
      ENDIF
      Jxp1 = Jx + 1
!
      CALL POS(Y,NY,Ywant,Jy)
      IF ( Jy==0 ) THEN
         Jy = 1
      ELSEIF ( Jy>=NY ) THEN
         Jy = NY - 1
      ENDIF
      Jyp1 = Jy + 1
!
      CALL POS(Z,NZ,Zwant,Jz)
      IF ( Jz==0 ) THEN
         Jz = 1
      ELSEIF ( Jz>=NZ ) THEN
         Jz = NZ - 1
      ENDIF
      Jzp1 = Jz + 1
!
      CALL POS(Q,NQ,Qwant,Jq)
      IF ( Jq==0 ) THEN
         Jq = 1
      ELSEIF ( Jq>=NQ ) THEN
         Jq = NQ - 1
      ENDIF
      Jqp1 = Jq + 1
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FINDINDEX4D
      SUBROUTINE FINDPOINT(Xwant,Ywant,Value,X,Y,Array,NX,NY,Jx,Jxp1,Jy,&
                         & Jyp1)
      IMPLICIT NONE
!
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1
      INTEGER :: NX , NY
      REAL :: Value , Xwant , Ywant
      REAL , DIMENSION(NX,NY) :: Array
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      INTENT (IN) Array , Jx , Jxp1 , Jy , Jyp1 , NX , NY , X , Xwant , &
                & Y , Ywant
      INTENT (OUT) Value
!
! Local variables
!
      REAL :: s1 , s2 , s3 , s4 , zt , zta , ztb , zu , zua , zub
!
!-----------------------------------------------------------------------
!
!     written by:   David Lary
!
!     started:      7/1/1993
!
!     last updated: 22/1/2004
!
!-----------------------------------------------------------------------
!
!     look up relevant grid points.
!
!-----------------------------------------------------------------------
!
      s1 = Array(Jx,Jy)
      s2 = Array(Jxp1,Jy)
      s3 = Array(Jxp1,Jyp1)
      s4 = Array(Jx,Jyp1)
!
!     find slopes used in interpolation;
!     i) x.
      zta = Xwant - X(Jx)
      ztb = X(Jxp1) - X(Jx)
!
      zt = zta/ztb
!
!     ii) y.
      zua = Ywant - Y(Jy)
      zub = Y(Jyp1) - Y(Jy)
!
      zu = zua/zub
!
!     use bilinear interpolation.
      Value = (1.0-zt)*(1.0-zu)*s1 + zt*(1.0-zu)*s2 + zt*zu*s3 +        &
            & (1.0-zt)*zu*s4
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FINDPOINT
      SUBROUTINE FINDPOINT2D(Xwant,Ywant,Value,X,Y,Array,NX,NY,Jx,Jxp1, &
                           & Jy,Jyp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1
      INTEGER :: NX , NY
      REAL :: Value , Xwant , Ywant
      REAL , DIMENSION(NX,NY) :: Array
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      INTENT (IN) Array , Jx , Jxp1 , Jy , Jyp1 , NX , NY , X , Xwant , &
                & Y , Ywant
      INTENT (OUT) Value
!
! Local variables
!
      REAL :: s1 , s2 , s3 , s4 , zt , zta , ztb , zu , zua , zub
!
!-----------------------------------------------------------------------
!
      s1 = Array(Jx,Jy)
      s2 = Array(Jxp1,Jy)
      s3 = Array(Jxp1,Jyp1)
      s4 = Array(Jx,Jyp1)
!
!     find slopes used in interpolation;
!     i) x.
      zta = Xwant - X(Jx)
      ztb = X(Jxp1) - X(Jx)
!
      zt = zta/ztb
!
!     ii) y.
      zua = Ywant - Y(Jy)
      zub = Y(Jyp1) - Y(Jy)
!
      zu = zua/zub
!
!     use bilinear interpolation.
      Value = (1.0-zt)*(1.0-zu)*s1 + zt*(1.0-zu)*s2 + zt*zu*s3 +        &
            & (1.0-zt)*zu*s4
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FINDPOINT2D
      SUBROUTINE FINDPOINT3D(Xwant,Ywant,Zwant,Value,X,Y,Z,Array,NX,NY, &
                           & NZ,Jx,Jxp1,Jy,Jyp1,Jz,Jzp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jx , Jxp1 , Jy , Jyp1 , Jz , Jzp1
      INTEGER :: NX , NY , NZ
      REAL :: Value , Xwant , Ywant , Zwant
      REAL , DIMENSION(NX,NY,NZ) :: Array
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) Array , Jx , Jxp1 , Jy , Jyp1 , Jz , Jzp1 , NX , NY , &
                & NZ , X , Xwant , Y , Ywant , Z , Zwant
      INTENT (INOUT) Value
!
! Local variables
!
      REAL :: f1 , f2 , f3 , s1 , s2 , s3 , s4 , s5 , s6 , s7 , s8
!
!----------------------------------------------------------------------
!
      f1 = (Xwant-X(Jx))/(X(Jxp1)-X(Jx))
      f2 = (Ywant-Y(Jy))/(Y(Jyp1)-Y(Jy))
      f3 = (Zwant-Z(Jz))/(Z(Jzp1)-Z(Jz))
!
!----------------------------------------------------------------------
!
      s1 = Array(Jx,Jy,Jz)
      s2 = Array(Jxp1,Jy,Jz)
      s3 = Array(Jx,Jyp1,Jz)
      s4 = Array(Jxp1,Jyp1,Jz)
      s5 = Array(Jx,Jy,Jzp1)
      s6 = Array(Jxp1,Jy,Jzp1)
      s7 = Array(Jx,Jyp1,Jzp1)
      s8 = Array(Jxp1,Jyp1,Jzp1)
!
!----------------------------------------------------------------------
!
      Value = 0.0
      Value = Value + (1.0-f1)*(1.0-f2)*(1.0-f3)*s1
      Value = Value + f1*(1.0-f2)*(1.0-f3)*s2
      Value = Value + (1.0-f1)*f2*(1.0-f3)*s3
      Value = Value + f1*f2*(1.0-f3)*s4
      Value = Value + (1.0-f1)*(1.0-f2)*f3*s5
      Value = Value + f1*(1.0-f2)*f3*s6
      Value = Value + (1.0-f1)*f2*f3*s7
      Value = Value + f1*f2*f3*s8
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FINDPOINT3D
      SUBROUTINE FINDPOINT4D(Xwant,Ywant,Zwant,Qwant,Value,X,Y,Z,Q,     &
                           & Array,NX,NY,NZ,NQ,Jx,Jxp1,Jy,Jyp1,Jz,Jzp1, &
                           & Jq,Jqp1)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Jq , Jqp1 , Jx , Jxp1 , Jy , Jyp1 , Jz , Jzp1
      INTEGER :: NQ , NX , NY , NZ
      REAL :: Qwant , Value , Xwant , Ywant , Zwant
      REAL , DIMENSION(NX,NY,NZ,NQ) :: Array
      REAL , DIMENSION(NQ) :: Q
      REAL , DIMENSION(NX) :: X
      REAL , DIMENSION(NY) :: Y
      REAL , DIMENSION(NZ) :: Z
      INTENT (IN) Array , Jq , Jqp1 , Jx , Jxp1 , Jy , Jyp1 , Jz ,      &
                & Jzp1 , NQ , NX , NY , NZ , Q , Qwant , X , Xwant , Y ,&
                & Ywant , Z , Zwant
      INTENT (INOUT) Value
!
! Local variables
!
      REAL :: f1 , f2 , f3 , f4 , s1 , s10 , s11 , s12 , s13 , s14 ,    &
            & s15 , s16 , s2 , s3 , s4 , s5 , s6 , s7 , s8 , s9
!
!----------------------------------------------------------------------
!
      f1 = (Xwant-X(Jx))/(X(Jxp1)-X(Jx))
      f2 = (Ywant-Y(Jy))/(Y(Jyp1)-Y(Jy))
      f3 = (Zwant-Z(Jz))/(Z(Jzp1)-Z(Jz))
      f4 = (Qwant-Q(Jq))/(Q(Jqp1)-Q(Jq))
!
!----------------------------------------------------------------------
!
      s1 = Array(Jx,Jy,Jz,Jq)
      s2 = Array(Jxp1,Jy,Jz,Jq)
      s3 = Array(Jx,Jyp1,Jz,Jq)
      s4 = Array(Jxp1,Jyp1,Jz,Jq)
      s5 = Array(Jx,Jy,Jzp1,Jq)
      s6 = Array(Jxp1,Jy,Jzp1,Jq)
      s7 = Array(Jx,Jyp1,Jzp1,Jq)
      s8 = Array(Jxp1,Jyp1,Jzp1,Jq)
      s9 = Array(Jx,Jy,Jz,Jqp1)
      s10 = Array(Jxp1,Jy,Jz,Jqp1)
      s11 = Array(Jx,Jyp1,Jz,Jqp1)
      s12 = Array(Jxp1,Jyp1,Jz,Jqp1)
      s13 = Array(Jx,Jy,Jzp1,Jqp1)
      s14 = Array(Jxp1,Jy,Jzp1,Jqp1)
      s15 = Array(Jx,Jyp1,Jzp1,Jqp1)
      s16 = Array(Jxp1,Jyp1,Jzp1,Jqp1)
!
!----------------------------------------------------------------------
!
      Value = 0.0
      Value = Value + (1.0-f1)*(1.0-f2)*(1.0-f3)*(1.0-f4)*s1
      Value = Value + f1*(1.0-f2)*(1.0-f3)*(1.0-f4)*s2
      Value = Value + (1.0-f1)*f2*(1.0-f3)*(1.0-f4)*s3
      Value = Value + f1*f2*(1.0-f3)*(1.0-f4)*s4
      Value = Value + (1.0-f1)*(1.0-f2)*f3*(1.0-f4)*s5
      Value = Value + f1*(1.0-f2)*f3*(1.0-f4)*s6
      Value = Value + (1.0-f1)*f2*f3*(1.0-f4)*s7
      Value = Value + f1*f2*f3*(1.0-f4)*s8
      Value = Value + (1.0-f1)*(1.0-f2)*(1.0-f3)*f4*s9
      Value = Value + f1*(1.0-f2)*(1.0-f3)*f4*s10
      Value = Value + (1.0-f1)*f2*(1.0-f3)*f4*s11
      Value = Value + f1*f2*(1.0-f3)*f4*s12
      Value = Value + (1.0-f1)*(1.0-f2)*f3*f4*s13
      Value = Value + f1*(1.0-f2)*f3*f4*s14
      Value = Value + (1.0-f1)*f2*f3*f4*s15
      Value = Value + f1*f2*f3*f4*s16
!
!----------------------------------------------------------------------
!
      END SUBROUTINE FINDPOINT4D
      SUBROUTINE POS(Xx,N,X,J)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: J
      INTEGER :: N
      REAL :: X
      REAL , DIMENSION(N) :: Xx
      INTENT (IN) N , X , Xx
      INTENT (OUT) J
!
! Local variables
!
      INTEGER :: jl , jm , ju
!
!-----------------------------------------------------------------------
!
!     Argrument list.
!
!     Name      Type              Description.
!     XX        Array of real     Monotonic array of length N.
!     (Unchanged on exit).
!
!     N         Integer           Length of array XX.
!     (Unchanged on exit).
!
!     X         Real              Value whose position in XX is
!     required.
!     (Unchanged on exit).
!
!     J         Integer           Index of X in array XX.
!     (Contains answer on exit).
!
!-----------------------------------------------------------------------
!
!     Given an array XX of length N, and given a value X, POS returns a
!     value J sucha that X lies between XX(J) and XX(J+1). XX must be
!     monotonic, either increasing or decreasing. J=0 or J=N is returned
!     to indicate that X is out of range.
!
!     The table entry J is found by bisection.
!
!     Based on Numerical Recipes, The art of scientific computing,
!     section 3.4, by Press, Flannery, Teukolsky & Vetterling,
!     Cambridge University Press, 1987.
!
!     Modified by: David Lary
!     ----------
!
!     Date Started : 7/2/1990
!
!     Last modified: 27/9/1991
!
!-----------------------------------------------------------------------
!
!     Initialize upper & lower limits.
      jl = 0
      ju = N + 1
!
!----------------------------------------------------------------------
!
      DO WHILE ( .TRUE. )
!
!----------------------------------------------------------------------
!
         IF ( .NOT..TRUE. ) THEN
            RETURN
!
         ELSEIF ( ju-jl>1 ) THEN
!
            jm = (ju+jl)/2
            IF ( Xx(N)>Xx(1) .EQV. X>Xx(jm) ) THEN
               jl = jm
            ELSE
               ju = jm
            ENDIF
!
!           Repeat untill the test condition is satisfied.
            CYCLE
         ENDIF
!
!        Set the output.
         J = jl
         EXIT
!
!----------------------------------------------------------------------
!
      ENDDO
!
!----------------------------------------------------------------------
!
      END SUBROUTINE POS

