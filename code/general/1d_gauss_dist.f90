
       program test_gauss

       implicit none
       integer i, nt
       real t(100), gauss(100)
       real mu, sig
       real, external :: gauss_dist_1d

       nt=100
       t(:) = 0.
       gauss(:) = 0.
       mu = 4.
       sig = 2.

       !call gauss_1d_dist(t,gauss,nt,mu,sig)
       gauss = gauss_dist_1d(t,nt,mu,sig)

       do i=1,nt
        print *, t(i), gauss(i)
       enddo


       endprogram test_gauss


       subroutine gauss_1d_dist(t,gauss,nt,mu,sig)

       !This subroutine computes a
       !normalized 1D gauss distribution

       !Imput:
       !nt:  length of the series
       !t:   serie to compute gaussian PDF
       !mu:  mean of the PDF
       !sig: sigma of the PDF

       !Output:
       !gauss: PDF

       implicit none
       integer i
       integer, intent(inout) :: nt
       real, intent(inout) :: t(nt), mu, sig
       real, intent(inout) :: gauss(nt)
       real dt, ti
       real pi/3.1415926535897932/
 
       ti = 0.
       dt = 0.1
       t(1) = ti+dt
       do i=2, nt
        t(i) = t(i-1) + dt
       enddo

       do i=1, nt
!         gauss_dist_1d(i) = 1./ (sig*sqrt(pi))* & 
!  &               exp(-1.*((real(i)-1.)*dt-mu)**2./sig**2.)
       enddo


       endsubroutine gauss_1d_dist



contains
       FUNCTION gauss_dist_1d(t,nt,mu,sig)


       implicit none
       integer i
       integer, intent(inout) :: nt
       real, intent(inout) :: t(nt), mu, sig
       real, dimension(nt) :: gauss_dist_1d(nt)
       real dt, ti
       real pi/3.1415926535897932/
 
       ti = 0.
       dt = 0.1
       t(1) = ti+dt
       do i=2, nt
        t(i) = t(i-1) + dt
       enddo

       do i=1, nt
         gauss_dist_1d(i) = 1./ (sig*sqrt(pi))* & 
  &               exp(-1.*((real(i)-1.)*dt-mu)**2./sig**2.)
       enddo



       END FUNCTION

