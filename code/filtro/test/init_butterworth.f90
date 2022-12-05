
!----------------------------------------------------------------------
!
!     Intialize the filter coefficients
!     The Butterworth filter was taken from:
!     http://jean-pierre.moreau.pagesperso-orange.fr/f_signal.html
!
!     IN : dt
!
!     IN/OUT : butt
!
!----------------------------------------------------------------------
MODULE INIT_BUTTERWORTH_MOD


CONTAINS

  SUBROUTINE INIT_BUTTERWORTH(butt, dt)

    !-----------------------     
    ! Subroutine parameters 
    !----------------------- 
    INCLUDE 'proc.h'
    TYPE(butter) ::  butt
    REAL dt

    !-----------------------------------
    ! Initialize the Butterworth filter
    !-----------------------------------
    ! Calculate the filter coefficients
    CALL Butterworth(butt%fc, dt, butt%order, butt%C, butt%NSections, butt%Tg)
    ! Initialize the filter coefficients
    CALL Init(0.0,butt%C,butt%NSections,butt%D)

    RETURN
  END SUBROUTINE INIT_BUTTERWORTH


!**************************************************************************
!*                       INIT FILTER PROCEDURE                            *
!* ---------------------------------------------------------------------- *
!* The filter response is initialized to stationnary value for a constant *
!* input signal value.                                                    *
!*                                                                        *
!* Calling mode:   INIT(Xdc,C,NSections,D);                               *
!* ---------------------------------------------------------------------- *
!* INPUTS:                                                                *
!* ------                                                                 *
!*        Xdc......: constant input value (double)                        *
!*        C........: Table[1..5,1..NSections] of filter coefficients      *
!*                   calculated previously by BUTTERWORTH subroutine      *
!*        NSections: Number of required 2nd order sections (integer)      *
!*                   = n/2     for n even (n=order of filter)             *
!*                   = (n+1)/2 for n odd                                  *
!*                   calculated previously by BUTTERWORTH subroutine      *
!* ---------------------------------------------------------------------- *
!* OUTPUTS:                                                               *
!* -------                                                                *
!*        D........: Table[1..2,1..NSections] of coefficients defining    *
!*                   the filter memory, initialized by INIT subroutine    *
!**************************************************************************
Subroutine Init(Xdc,C,NSections,D)
real Xdc,C(5,10), D(2,10)
real dc,Csum 
    dc=Xdc
    DO j=1, NSections
      D(2,j)=dc/(1.d0-C(1,j)-C(2,j))
      D(1,j)=D(2,j)
      Csum=0.d0
      do ii=1, 4
        Csum=Csum + C(ii,j)
      end do
      dc=C(5,j)*(dc+D(2,j)*Csum)
    END DO
    RETURN
End subroutine

!**********************************************************************
!*          Calculates the Butterworth filter coefficients            *
!* ------------------------------------------------------------------ *
!*  Calling mode:   Butterworth(Fc,Ts,n,C,NSections,Tg);              *
!* ------------------------------------------------------------------ *
!*  INPUTS:                                                           *
!*  ------                                                            *
!*         Fc.......: Cut off frequency                               *
!*         dt.......: Sampling time of input signal                   *
!*         N........: Order of filter (1 to 4)                        *
!* ------------------------------------------------------------------ *
!*  OUTPUTS:                                                          *
!*  -------                                                           *
!*         C........: Table[1..5,1..NSections] of filter coefficients *
!*                    calculated previously by BUTTERWORTH subroutine *
!*         NSections: Number of required 2nd order sections (integer) *
!*                    = N/2     for N even                            *
!*                    = (N+1)/2 for N odd                             *
!*         Tg.......: Group delay in seconds                          *
!*********************************************************************}
Subroutine Butterworth(Fc,dt,N,C,NSections,Tg)
real Fc,dt,C(5,10),Tg
real arg,m,Omega,OmegaSq,Pi,temp,xmodn,xn,xns2

  Zero = 0.d0
  ONE  = 1.d0
  TWO  = 2.d0
  HALF = 0.5d0

  Pi = 3.1415926535d0
  arg=Pi*dt*Fc
!  If (dabs(arg) > 2.d0*Pi) THEN
  If (abs(arg) > 2.d0*Pi) THEN
    m=INT(arg/2.d0/Pi)
    arg=arg-(m*2.d0*Pi)
  END IF
!  Omega= dtan(arg)
  Omega= tan(arg)
  OmegaSq=Omega*Omega
  xn=N
  If (xn/2.eq.INT(xn/2)) THEN
    xmodn=Zero; temp=HALF
  else
    xmodn=HALF; temp=Zero
  end if
!HUGO CHANGED N/2 TO N/2. TO AVOID ERRORS WHEN N=1
  xns2 = float(N)/2.
  ns2=N/2
  NSections=INT(xns2+xmodn)
  Tg=Zero
  If (N>1) THEN
    DO i=1, ns2
!      Rep=Omega*DCOS(Pi*(i-temp)/N)
      Rep=Omega*COS(Pi*(i-temp)/N)
      Tg=Tg+dt*Rep/OmegaSq
      W0=TWO*Rep
      W1=ONE + W0+OmegaSq
      C(1,i)=-TWO*(OmegaSq-ONE)/W1
      C(2,i)=-(ONE-W0+OmegaSq)/W1
      C(3,i)=TWO
      C(4,i)=ONE
      C(5,i)=OmegaSq/W1
    end do
  end if
  If (temp.eq.Zero) THEN
    C(1,Nsections)=(ONE-Omega)/(ONE+Omega)
    C(2,NSections)= Zero
    C(3,NSections)= ONE
    C(4,NSections)= Zero
    C(5,NSections)= Omega/(ONE+Omega)
    Tg= Tg+dt/(TWO*Omega)
  END IF
RETURN 
End subroutine !of Butterworth()

END MODULE INIT_BUTTERWORTH_MOD

