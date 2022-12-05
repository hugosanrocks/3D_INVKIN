      program probing

      implicit none
      !Counters
      integer :: i, j, k, l, n_samp, iunit
      integer :: iunit1, iunit2, iunit3
      !Size of arrays
      integer :: mm, nn, m, n, cont, cont2
      !Size of fault plane
      integer :: n_stk, n_dip, n_time
      !Coefficients for matrix multiplication
      real :: alpha, beta, mean, sigma
      !Arrays for matrix and vectors
      real, dimension(:,:), allocatable :: hess, vec, xsol
      !parameters
      integer :: corrsize, tsamw, n_win, sam_res, realiz,&
  &              corrsizex, corrsizez, xsamw, zsamw, n_winx,&
  &              n_winz, sam_resx, sam_resz
      real :: dt, dx, dz, twin, tfin, res, xwin, zwin
      !correlation arrays
      real, dimension(:), allocatable :: tserie1, tserie2, tserie3
      real, dimension(:), allocatable :: array_x, array_z
      real, dimension(:), allocatable :: x_1, x_2, x_corr
      real, dimension(:), allocatable :: z_1, z_2, z_corr


      !Number of nodes
      n_stk = 24
      n_dip = 12
      n_time = 57

      !Number of realizations
      realiz = 50

      !dt of signals
      dt = 0.16
      tfin = 8.96
      twin = 2.
      !Dimensions along space
      dx = 1.5
      dz = 1.5
      !length of window along space
      xwin = 5.
      zwin = 5.

      !Size of Hessian matrix MxM
      m = n_stk*n_dip*n_time  !siv1=288*57 !toysmall4x4 100*16
      !Size of random model vector Mx1
      n = m

      if ( m .ne. m ) then
        write(*,*) ' **ERROR M .ne. N** '
      endif

      !Allocate memory for Hessian and vector
      allocate(hess(m,m),vec(n,1),xsol(n,1))
      allocate(array_x(n),array_z(n))

      !Read Hessian information
      iunit = 22
      !If compiled with ifort recl=m
      !If compiled with gfortran recl=m*4
      open(iunit,file='hess.full',&
  &        status='old',action='read',&
  &        form='unformatted',access='direct',&
  &        recl=m)
      do i=1,m
        read(iunit,rec=i) hess(:,i)
      enddo
      close(iunit)

      !Prepare random vector with Gaussian dist
      mean = 0.
      sigma = 0.1

      iunit=25
      open(iunit,file='Output_data.out',status='unknown')
      write(iunit,*) mean, sigma, ' Mean and Sigma used'

      !Matrix multiplication constants
      nn = 1
      alpha = 1.
      beta = 0.

      !Samples per time window
      tsamw = int(twin/dt) + 1  !+1 is the t=0
      !Size of correlation array
      corrsize = tsamw + tsamw - 1
      !Number of time windows per node
      n_win = int(n_time/tsamw)
      !Number of samples out of the window
      sam_res = n_time - (n_win*tsamw)

      !Space correlation windows
      !Samples per space window
      xsamw = int(xwin/dx) + 1  !+1 is the x=0
      zsamw = int(zwin/dz) + 1  !+1 is the x=0
      !Size of correlation array
      corrsizex = (xsamw*2) - 1
      corrsizez = (zsamw*2) - 1
      !Number of time windows per node
      n_winx = int(n_stk/xsamw)
      n_winz = int(n_dip/zsamw)
      !Number of samples out of the window
      sam_resx = n_stk - (n_winx*xsamw)
      sam_resz = n_dip - (n_winz*zsamw)

      !Print values only to check
      write(25,*) tsamw, n_win, sam_res
      write(25,*) 'Time samp wind, no. wind per node, samp res'
      write(25,*) xsamw, n_winx, sam_resx
      write(25,*) 'Stk samp wind, no. wind per node, samp res'
      write(25,*) zsamw, n_winz, sam_resz
      write(25,*) 'Dip samp wind, no. wind per node, samp res'
      close(iunit)

      !Allocate memory for correlations
      allocate(tserie1(tsamw),tserie2(tsamw),tserie3(corrsize))
      allocate(x_1(xsamw),x_2(xsamw),x_corr(corrsizex))
      allocate(z_1(xsamw),z_2(xsamw),z_corr(corrsizex))


      iunit1=21
      open(iunit1,file='Corr_time.out',status='unknown')
      iunit2=22
      open(iunit2,file='Corr_strike.out',status='unknown')
      iunit3=23
      open(iunit3,file='Corr_dip.out',status='unknown')



   do n_samp=1,realiz
      vec(:,:) = 0.

      !Create random model vector
      call random_vector(mean,sigma,vec,m)

      !Hessian vector product
      call sgemm('N','N',m,nn,m,alpha,hess,m,vec,m,beta,xsol,m)

      !Flush arrays
      array_x(:) = 0.
      array_z(:) = 0.
      !Arrange vector according to strike and dip
      l = 1
      do i=1,n_time
       cont = i
       do j=1,n_dip
        do k=1,n_stk
         cont2=(k-1)*n_dip+(j-1)+((i-1)*(n_dip*n_stk))+1
         !strike
         array_x(l) = xsol(cont,1)
         !dip
         array_z(cont2) = xsol(cont,1)
         cont = cont + n_time
         l = l + 1
        enddo
       enddo
      enddo

      !In case dip arragement is wrong
      !Arrange vector according to dip
      !l = 1
      !do i=1,2!n_time
      ! do j=1,n_stk
      !  do k=1,n_dip
      !   cont=(j-1)*57+57*24*(k-1)+i
      !   array_z(l) = xsol(cont,1)
      !   l = l + 1
      !  enddo
      ! enddo
      !enddo

      !Correlations in time
      !Extract each window from the H*v
      k = 1
      do i=1,n_stk*n_dip
        do j=1,n_win
          tserie1(1:tsamw) = xsol(k:k+tsamw-1,1)  !Extracting the window
          tserie2(:) = tserie1(:)
          tserie3(:) = 0.
          call xcorr1d(tserie1,tsamw,tserie2,tsamw,&
  &                    tserie3,corrsize)
          do l=1,corrsize
           write(iunit1,*) tserie3(l)
          enddo
          k = k+tsamw    !jump to next window
        enddo
        k = k+sam_res    !jump ramaining samples
      enddo

      !Correlations along strike
      !Extract each window from the H*v
      k = 1
      do i=1,n_dip*n_time
        do j=1,n_winx
          x_1(1:xsamw) = array_x(k:k+xsamw-1) !Extracting the window
          x_2(:) = x_1(:)
          x_corr(:) = 0.
          call xcorr1d(x_1,xsamw,x_2,xsamw,&
  &                    x_corr,corrsizex)
          do l=1,corrsizex
           write(iunit2,*) x_corr(l)
          enddo
          k = k+xsamw    !jump to next window
        enddo
        k = k+sam_resx    !jump ramaining samples
      enddo

      !Correlations along dip
      !Extract each window from the H*v
      k = 1
      do i=1,n_stk*n_time
        do j=1,n_winz
          z_1(1:zsamw) = array_z(k:k+zsamw-1)  !Extracting the window
          z_2(:) = z_1(:)
          z_corr(:) = 0.
          call xcorr1d(z_1,zsamw,z_2,zsamw,&
  &                    z_corr,corrsizez)
          do l=1,corrsizez
           write(iunit3,*) z_corr(l)
          enddo
          k = k+zsamw    !jump to next window
        enddo
        k = k+sam_resz    !jump ramaining samples
      enddo

      write(*,*) ' Finished sample no. ', n_samp

   enddo
   close(iunit1)
   close(iunit2)
   close(iunit3)

      !write to check
      do i=1,m
       write(23,*) vec(i,1), xsol(i,1)
       write(26,*) hess(i,7430)
      enddo

      !Deallocate memory
      deallocate(hess,vec,xsol)
      deallocate(tserie1,tserie2,tserie3)
      deallocate(array_x,array_z)
      deallocate(x_1,x_2,z_1,z_2)
      endprogram probing
