
      subroutine tractfft (proc_mesh,tract,n,nc,nrec)

      implicit none
      integer, intent(inout) :: n, nc, nrec
      real, intent(inout) :: tract(3,n)

      !Define all the variables needed to read stresses and calculate the tractions
      !associated "Green's functions", variables in include/green.h
      INCLUDE 'proc.h'
      TYPE (mesh) :: proc_mesh
      include "fftw3.f"

      integer j, k, cont
      integer*8 plan_forward
!      double precision in(n), in2(n)
!      double complex out(nc)
      real in(n)
      complex out(nc)
      integer iunit, reclent


      INQUIRE(iolength=reclent) out(1)

      call sfftw_plan_dft_r2c_1d_ ( plan_forward, n, in, out,&
     &  FFTW_ESTIMATE )

      iunit=16
      OPEN(iunit,file='dat/TRACT_fft.bin',status='unknown',&
 &         form='unformatted',ACCESS='DIRECT',recl=reclent*nc)

      cont=1 + ((nrec-1)*proc_mesh%ncomp)
        in(:)=0.
        do k = 1, proc_mesh%ncomp
        do j = 1, proc_mesh%interp_i
          in(j) = tract(k,j)
        enddo

       ! do i = 1, n
       !   write ( 4, '(2x,i4,2x,g14.6)' ) i, in(i)
       ! end do

!  Set up a plan, and execute the plan to transform the IN data to
!  the OUT FFT coefficients.

        call sfftw_execute_ ( plan_forward )

       ! do i = 1, nc
       !   write ( 4, '(2x,i4,2x,2g14.6)' ) i, out(i)
       ! end do

        write(iunit,rec=cont) out(:)
        print *, cont, proc_mesh%sta_i, proc_mesh%comp_i, 'freq'

        !print *, 'rec=saved',cont
     ! call sfftw_plan_dft_c2r_1d_ ( plan_backward, n, out, in2,&
     !&  FFTW_ESTIMATE )

     ! call sfftw_execute_ ( plan_backward )

     ! write ( 4, '(a)' ) ' '
     ! write ( 4, '(a)' ) '  Recovered input data divide by N:'
     ! write ( 4, '(a)' ) ' '

     ! do i = 1, n
     !   write ( 4, '(2x,i4,2x,g14.6)' ) i, in2(i) / dble ( n )
     ! end do

          cont = cont + 1
      enddo

      call sfftw_destroy_plan_ ( plan_forward )
     ! call dfftw_destroy_plan_ ( plan_backward )
      close(iunit)



      return
      end subroutine tractfft
