c     Simple program to demonstrate calling the wrapper routines
c     to perform 1D transforms in Fortran.  This program should be
c       linked with -lfftw -lm.

      program test

      implicit none

      include "fftw3.f"

      integer N
      parameter(N=4)

      integer*8 plan
      double complex in, out
      dimension in(N),out(N)

      integer i

      write(*,*) 'Input array:'

      do i = 1,N,1
        in(i) = dcmplx(float(i),float(i+1))
        write(*,*) '    in(',i,') = ',in(i)
      enddo

      call dfftw_plan_dft_1d ( plan, N, in, out,
     ^                         FFTW_FORWARD, FFTW_ESTIMATE )

      call dfftw_execute ( plan )

      write(*,*) 'Output array:'
      do i = 1,N,1
        write(*,*) '    out(',i,') = ',out(i)
      enddo

      call dfftw_destroy_plan ( plan )

      call dfftw_plan_dft_1d ( plan, N, out, in,
     ^                         FFTW_FORWARD, FFTW_ESTIMATE )

      call dfftw_execute ( plan )

      write(*,*) 'Output array after inverse FFT:'
      do i = 1,N,1
        write(*,*) '    ',N,' * in(',i,') = ',in(i)
      enddo

      call dfftw_destroy_plan ( plan )

      end
