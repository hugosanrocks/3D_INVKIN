!----------------------------------------------------------------------

!----------------------------------------------------------------------

MODULE INIT_RANDOM_SEED_MOD
 
  ! This module is only required when using Intel compiler
USE IFPORT

CONTAINS

  SUBROUTINE INIT_RANDOM_SEED()
    !----------------------------------------------------------------------
    !
    !  INIT_RANDOM_SEED
    !
    !  Initialize the random seed with a varying seed in order to ensure 
    !  a different random number sequence for each invocation of the program
    ! 
    !  Taken from: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    !----------------------------------------------------------------------   

    IMPLICIT NONE 
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
    INTEGER(8) :: count, tms
          
    CALL random_seed(size = n)
    ALLOCATE(seed(n))
    ! First try if the OS provides a random number generator
    OPEN(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    IF (istat == 0) THEN
       READ(un) seed
       CLOSE(un)
    ELSE
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       CALL system_clock(count)
       IF (count /= 0) THEN
          t = transfer(count, t)
       ELSE
          CALL date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       END IF
       s = ieor(t(1), t(2))
       pid = GETPID() + 1099279 ! Add a prime
       s = ieor(s, pid)
       IF (n >= 3) THEN
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          IF (n > 3) THEN
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          END IF
       ELSE
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       END IF
    END IF
    CALL random_seed(put=seed)

    RETURN
  END SUBROUTINE INIT_RANDOM_SEED
END MODULE INIT_RANDOM_SEED_MOD
