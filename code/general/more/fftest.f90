        program test

        implicit none
       
        call test01( )

        end program test



      subroutine test01

!*********************************************************************72
!
!c MAIN is the main program for FFTW3_PRB.
!
!  Discussion:
!
!    FFTW3_PRB tests the FFTW3 library.
!
!  Modified:
!
!    05 November 2007
!
!  Author:
!
!    John Burkardt
!
!*********************************************************************72
!
!c TEST02: real 1D data.
!
!  Modified:
!
!    23 October 2005
!
      implicit none

      integer n
      integer nc
      integer nconv

      parameter ( n = 10 )
      parameter ( nc = 11 )
      parameter ( nconv = 20 )

      include "fftw3.f"

!      double precision r8_uniform_01
      integer i
      real in(nconv), iny(nconv), inh(nconv)
      real in2(nconv), in2y(nconv)
      complex out(nc), outy(nc), hf(nc)
      integer*8 plan_backward
      integer*8 plan_forward
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFT1D_Real'


      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Demonstrate FFTW3 on a single vector'
      write ( *, '(a)' ) '  of real data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Transform data to FFT coefficients.'
      write ( *, '(a)' ) '  Backtransform FFT coefficients to recover '
      write ( *, '(a)' ) '  data.'
      write ( *, '(a)' ) '  compare recovered data to original data.'
!
!  Set up the input data, a real vector of length N.
!
!      do i = 1, n
!        in(i) = r8_uniform_01 ( seed )
!      end do

         in(1) = 1.
         in(2) = 2.
         in(3) = 4.
         in(4) = 5.
         in(5) = 2.
         in(6) = 3.
         in(7) = 2.
         in(8) = 1.
         in(9) = 0.
         in(10) = 3.
         do i = 11 , 20
           in(i) = 0.
         enddo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, in(i)
      end do


!
!  Set up a plan, and execute the plan to transform the IN data to
!  the OUT FFT coefficients.
!
      call sfftw_plan_dft_r2c_1d_ ( plan_forward, nconv, in, out,&
     &  FFTW_ESTIMATE )

      call sfftw_execute_ ( plan_forward )


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nc
        write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
      end do



!
!  Set up a plan, and execute the plan to backtransform the
!  complex FFT coefficients in OUT to real data.


      call sfftw_plan_dft_c2r_1d_ ( plan_backward, nconv, out, inh,&
     &  FFTW_ESTIMATE )

      call sfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, inh(i) / dble ( nconv )
      end do

!
!  Release the memory associated with the plans.
!
      call sfftw_destroy_plan_ ( plan_forward )
      call sfftw_destroy_plan_ ( plan_backward )

!
!  Terminate.
!
      return
      end


