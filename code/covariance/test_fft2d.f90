       program test_fft2

        implicit none
        include "fftw3.f"
        double precision in
        dimension in(4,3)
        double complex out
        dimension out(3, 3)
        integer*8 plan, i, j, m, n

        m=4
        n=3

        do i=1,4
          do j=1,3
            in(i,j) = real(i+j)
          enddo
          print *, in(i,:)
        enddo

        call dfftw_plan_dft_r2c_2d(plan,M,N,in,out,FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan, in, out)
        call dfftw_destroy_plan(plan)

        do i=1,3
          print *,out(i,:)
        enddo




       endprogram test_fft2
