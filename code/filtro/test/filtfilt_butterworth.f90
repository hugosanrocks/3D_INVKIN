!----------------------------------------------------------------------
!
!     Apply a zero phase filtering with a Butterworth filter 
!     The Butterworth filter was taken from:
!     http://jean-pierre.moreau.pagesperso-orange.fr/f_signal.html
!
!     IN : butt, nstep
!
!     IN/OUT : tseries
!
!----------------------------------------------------------------------
MODULE FILTFILT_BUTTERWORTH_MOD

CONTAINS

  SUBROUTINE FILTFILT_BUTTERWORTH(tseries, butt, nstep)

    !-----------------------     
    ! Subroutine parameters 
    !----------------------- 
    REAL, DIMENSION(:), POINTER :: tseries
    INCLUDE 'proc.h'
    TYPE(butter) ::  butt
    INTEGER nstep

    !-----------------------
    ! Local variables
    !-----------------------
    REAL, DIMENSION(:), POINTER :: ftseries
    INTEGER ii   

    ! Allocate the required memory 
    ALLOCATE(ftseries(nstep))

    !-------------------
    ! Forward filtering
    !-------------------    
    ! Recursively call Butterworth filter
    DO ii = 1,nstep
       CALL Filter(ftseries(ii),tseries(ii),butt%NSections,butt%C,butt%D)
    ENDDO

    ! Reverse filtered time series
    DO ii = 1,nstep
       tseries(ii) = ftseries(nstep+1-ii)
    ENDDO

    !-------------------
    ! Reverse filtering
    !-------------------  
    ! Apply again the filter
    DO ii = 1,nstep
       CALL Filter(ftseries(ii),tseries(ii),butt%NSections,butt%C,butt%D)
    ENDDO

    ! Reverse filtered time series
    DO ii = 1,nstep
       tseries(ii) = ftseries(nstep+1-ii)
    ENDDO  

    ! Deallocate the requested memory
    DEALLOCATE(ftseries)   

    RETURN
  END SUBROUTINE FILTFILT_BUTTERWORTH

!***********************************************************************
!*          Filtering a signal F(t) by Butterworth method              *
!*             (removing frequencies greater then Fc)                  *
!* ------------------------------------------------------------------- *
!* Calling mode:   Filter(Xs,Xd,NSections,C,D);                        *
!* ------------------------------------------------------------------- *
!* INPUTS:                                                             *
!* -------                                                             *
!*      Xd.......: current value of input signal (double)              *
!*      NSections: Number of required 2nd order sections (integer)     *
!*                   = n/2     for n even                              *
!*                   = (n+1)/2 for n odd                               *
!*      C........: Table[1..5,1..NSections] of filter coefficients     *
!*                 calculated previously by BUTTERWORTH subroutine     *
!*      D........: Table[1..2,1..NSections] of coefficients defining   *
!*                 the filter memory, initialized by INIT subroutine   *
!* ------------------------------------------------------------------- *
!* OUTPUTS:                                                            *
!* -------                                                             *
!*      D..........: Table updated after the call to Filter subroutine *
!*      Xs.........: current value of filtered signal (double)         *
!* ------------------------------------------------------------------- *
!* Reference                                                           *
!* ---------                                                           *
!*  "Lawrence R.Rabiner et Bernard Gold                                *
!*   Theory and application of digital processing.                     *
!*   Prentice Hall Inc., EnglewoodclIFfs,NEW JERSEY,1975."             *
!*   [BIBLI 15].                                                       *
!*                                                                     *
!*                                   F90 version by J-P Moreau, Paris  *
!*               from Fortran version by J-P Dumont / Dang Trong Tuan  *
!***********************************************************************
Subroutine Filter(Xs,Xd,NSections,C,D)
real  Xs,Xd,C(5,10), D(2,10)
    x=Xd
    do ii=1, NSections
      xerr=x+C(1,ii)*D(1,ii)+C(2,ii)*D(2,ii)
      y=C(5,ii)*(xerr +C(3,ii)*D(1,ii)+C(4,ii)*D(2,ii))
      D(2,ii)=D(1,ii)
      D(1,ii)=xerr
      x=y
    end do
    Xs=x
    return
End subroutine

END MODULE FILTFILT_BUTTERWORTH_MOD

