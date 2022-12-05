! --------------------------------------------------------------------
!       SUBROUTINE TO GENERATE PSEUDORANDOM NORMAL NUMBERS
! --------------------------------------------------------------------

!----------------------------------------------------------------------
!
!  NORMAL_DIST
!
!  The normal probability distribution function (PDF) is sampled,
!  with mean A and standard deviation B.
! 
!  This subroutine is a modification of the function r8_normal_ab
!  Taken from: http://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.html
!----------------------------------------------------------------------   

SUBROUTINE NORMAL_DIST (a, b, x)

  IMPLICIT NONE

  ! Parameters
  REAL a ! Mean
  REAL b ! Standard Deviaton
  REAL x
 
  ! Local variables
  REAL, PARAMETER :: PI = 3.141592653589793D+00
  REAL r1, r2

  CALL RANDOM_NUMBER (HARVEST = r1)
  CALL RANDOM_NUMBER (HARVEST = r2)

  x = a + b * (SQRT ( - 2.0D+00 * LOG ( r1 ) ) * COS ( 2.0D+00 * PI * r2 ))

END SUBROUTINE NORMAL_DIST
