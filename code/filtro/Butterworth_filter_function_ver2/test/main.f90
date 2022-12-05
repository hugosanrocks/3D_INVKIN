!----------------------------------------------------------------------
!
!     Test of the Butterworth filter to be incorporated in the
!               adjoint kinematic inversion procedure
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
#include "precis.h"
!----------------------------------------------------------------------
PROGRAM MAIN

  USE INIT_BUTTERWORTH_MOD
  USE FILTFILT_BUTTERWORTH_MOD

  !----------------------------------------------------------------------
  INCLUDE 'param.h'
  !----------------------------------------------------------------------
  INCLUDE 'paral.h'
  !----------------------------------------------------------------------
#include "global.h"
  !----------------------------------------------------------------------
#include "type_def.h"
  !----------------------------------------------------------------------

  !------------------------------     
  !     variables declaration
  !------------------------------
  TYPE(butterworth_filter_type) butt 
  MYFLOAT, DIMENSION(:), POINTER :: tseries
  MYFLOAT, DIMENSION(:,:), POINTER :: copy_tseries
  MYFLOAT, DIMENSION(:), POINTER :: tseries_total
  INTEGER nstep
  MYFLOAT dt
  INTEGER ii, rec_sel, nrec
  ! For printing purposes
  MYFLOAT t, tend, tbegin

  ! Number of data in the seismograms
  nstep = 2188
  ! Number of receivers
  nrec = 80
  ! Receivers selection
  rec_sel = 1
  ! Time step
  dt = 0.005486
  ! Final and initial times (just for printing)
  tend = 12.002882
  tbegin = 0.0

  ! Allocate the required memory 
  ALLOCATE(tseries_total(nstep*nrec))
  ALLOCATE(tseries(nstep))
  ALLOCATE(copy_tseries(nstep,2))

  ! Filter parameters
  butt%order = 4
  butt%fc = 3.  

  ! Open the seismograms file
  OPEN(UNIT = 1, FILE = 'VZ', ACCESS='DIRECT', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', RECL = nstep*nrec*4)
  ! Read all the seismograms
  READ(1,rec=1) tseries_total
  ! Close the file
  CLOSE (UNIT=1)

  ! Initialize the butterworth filter
  CALL INIT_BUTTERWORTH(butt, dt)

  ! Read the seismogram of the selected receiver
  DO ii = 1,nstep
     tseries(ii) = tseries_total((ii-1)*nrec+1)
  ENDDO

  ! Apply the butterworth filter
  CALL FILTFILT_BUTTERWORTH(tseries, butt, nstep)
  copy_tseries(:,2) = tseries(:)

  ! Print the resultsprint input and filtered signals to output file
  OPEN(UNIT=1,FILE='test_butterworth.txt',STATUS='UNKNOWN')
  t=tbegin
  DO ii = 1,nstep
    WRITE(1,'(E12.5,x,E12.5,x,E12.5)') t, tseries_total((ii-1)*nrec+1), copy_tseries(ii,2)
    t=t+dt
  ENDDO
  CLOSE(unit=1)
  PRINT *,' '
  PRINT *,'Results in file test_butterworth.txt :)'

  ! Deallocate the required memory
  DEALLOCATE(tseries_total)
  DEALLOCATE(tseries)
  DEALLOCATE(copy_tseries)

  STOP
END PROGRAM MAIN
