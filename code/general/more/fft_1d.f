      subroutine test01

c*********************************************************************72
c
cc MAIN is the main program for FFTW3_PRB.
c
c  Discussion:
c
c    FFTW3_PRB tests the FFTW3 library.
c
c  Modified:
c
c    05 November 2007
c
c  Author:
c
c    John Burkardt
c
c*********************************************************************72
c
cc TEST02: real 1D data.
c
c  Modified:
c
c    23 October 2005
c
      implicit none

      integer n
      integer nc
      integer nconv

      parameter ( n = 10 )
      parameter ( nc = 11 )
      parameter ( nconv = 20 )

      include "fftw3.f"

c      double precision r8_uniform_01
      integer i
      double precision in(nconv), iny(nconv), inh(nconv)
      double precision in2(nconv), in2y(nconv)
      double complex out(nc), outy(nc), hf(nc)
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
      write ( *, '(a)' ) '  Compare recovered data to original data.'
c
c  Set up the input data, a real vector of length N.
c
c      do i = 1, n
c        in(i) = r8_uniform_01 ( seed )
c      end do

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

      do i= 1 , n
        iny(i) = in(i) * 2.
      end do
        iny(11) = 4.
        do i = 12 , 20
           in(i) = 0.
        enddo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, in(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Y Data:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, iny(i)
      end do

c
c  Set up a plan, and execute the plan to transform the IN data to
c  the OUT FFT coefficients.
c
      call dfftw_plan_dft_r2c_1d_ ( plan_forward, nconv, in, out, 
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nc
        write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
      end do



      call dfftw_plan_dft_r2c_1d_ ( plan_forward, nconv, iny, outy,
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )
      

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Y Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nc
        write ( *, '(2x,i4,2x,2g14.6)' ) i, outy(i)
      end do



c     CONVOLUTION

      do i = 1, nc
        hf(i) = out(i) * outy(i)
      enddo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Multiplication result:'
      write ( *, '(a)' ) ' '

      do i = 1, nc
        write ( *, '(2x,i4,2x,2g14.6)' ) i, hf(i)
      end do



c
c  Set up a plan, and execute the plan to backtransform the
c  complex FFT coefficients in OUT to real data.
c
      call dfftw_plan_dft_c2r_1d_ ( plan_backward, nconv, out, in2,
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data divide by N:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, in2(i) / dble ( nconv )
      end do




      call dfftw_plan_dft_c2r_1d_ ( plan_backward, nconv, outy, in2y,
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input y data divide by N:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, in2y(i) / dble ( nconv )
      end do






      call dfftw_plan_dft_c2r_1d_ ( plan_backward, nconv, hf, inh,
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered convolution result:'
      write ( *, '(a)' ) ' '

      do i = 1, nconv
        write ( *, '(2x,i4,2x,g14.6)' ) i, inh(i) / dble ( nconv )
      end do

c
c  Release the memory associated with the plans.
c
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

c
c  Terminate.
c
      return
      end

      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit double precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    23 October 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end


      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
