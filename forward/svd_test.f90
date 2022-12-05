      PROGRAM svd_test

      INTEGER          M, N
      PARAMETER        ( M = 6, N = 4 )
      INTEGER          LDA, LDU, LDVT
      PARAMETER        ( LDA = M, LDU = M, LDVT = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )

!     .. Local Scalars ..
      INTEGER          INFO, LWORK, i
!
!     .. Local Arrays ..
!     IWORK dimension should be at least 8*MIN(M,N)
      INTEGER          IWORK( 8*N )
      DOUBLE PRECISION A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ),&
  &                    WORK( LWMAX )
!      DATA             A/
!  &     7.52,-0.76, 5.13,-4.75, 1.33,-2.40,&
!  &    -1.10, 0.62, 6.62, 8.52, 4.91,-6.77,&
!  &    -7.95, 9.34,-5.66, 5.75,-5.49, 2.34,&
!  &     1.08,-7.10, 0.87, 5.30,-3.52, 3.95 &
!  &                     /
!
!     .. External Subroutines ..
!      EXTERNAL         DGESDD
!      EXTERNAL         PRINT_MATRIX
!
!     .. Intrinsic Functions ..
!      INTRINSIC        INT, MIN
!
!     .. Executable Statements ..
      WRITE(*,*)'DGESDD Example Program Results'
!
!     Query the optimal workspace.
!
      A(1:6,1) = [ 7.52,-0.76, 5.13,-4.75, 1.33,-2.40]
      A(1:6,2) = [-1.10, 0.62, 6.62, 8.52, 4.91,-6.77]
      A(1:6,3) = [-7.95, 9.34,-5.66, 5.75,-5.49, 2.34]
      A(1:6,4) = [ 1.08,-7.10, 0.87, 5.30,-3.52, 3.95]


      LWORK = -1
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
  &                LDVT, WORK, LWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!     Compute SVD.
!
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
  &                LDVT, WORK, LWORK, IWORK, INFO )
!
!     Check for convergence.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF
!
!     Print singular values.
!
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )

!     Print left singular vectors.
!
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',&
  &                      M, N, U, LDU )
!
!     Print right singular vectors.
!
!print *, U(1,:)
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',&
  &                      N, N, VT, LDVT )
      STOP
      END
!
!     End of DGESDD Example.
!
!  =============================================================================
!
!     Auxiliary routine: printing a matrix.
!
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END
