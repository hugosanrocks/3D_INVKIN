     subroutine identify_past(green_mesh)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      integer :: iunit, i, j, k, ii, jj, cont
      integer :: predsize
      real :: timewin, timemark, percent, newweight
      integer :: tail, fixw
      real,dimension(:),allocatable :: serie
      real,dimension(:,:),allocatable :: weight
      integer,dimension(:),allocatable :: idsub

      !Percetange of time history already explained
      percent = green_mesh%percent_win !0.7before change  !0.7-real(green_mesh%interp_i-13)/140.
      print *, ' percentage of frozen time window per node:', percent

      !New weight to freeze prior model to previous solutions
      newweight = green_mesh%new_weight
      !newweight = 4.before change
      print *, ' new weight used for frozen time percentage:', newweight

      allocate(serie(green_mesh%msub))
      if ( green_mesh%rake_opt .eq. 1 ) then
       allocate(idsub(green_mesh%msub))
       allocate(weight(green_mesh%interp_i,green_mesh%msub))
      else
       allocate(idsub(green_mesh%msub))
       allocate(weight(green_mesh%interp_i,green_mesh%msub))
      endif

      weight(:,:) = 0.
      idsub(:) = 0.

      iunit=25
      open(iunit,file='dat/times_rup.dat',status='unknown')
      do i=1,green_mesh%msub
       read(iunit,*) serie(i)
      enddo
      close(iunit)

      green_mesh%idsub(:) = 0              !Order of subfaults to predict

      !Time at the source location at every stage
      timewin = (real(green_mesh%samp4win(green_mesh%synwin))-1.)&
  &             *green_mesh%slipdt
      !Do not update anything before this mark
      timemark = timewin - percent*timewin
      print *, ' Fixed time not to update:', timemark

      !Sample where the tail starts
      tail = int(timemark/green_mesh%slipdt) + 2 !0 sample + 1 more
      print *, ' Tail to remove starts at: ', tail

      !Arrange weight in a matrix similar to the slip-rate
      k=1
      do i=1,green_mesh%interp_i
       do j=1,green_mesh%msub
        weight(i,j) = green_mesh%diag(k)
        k=k+1
       enddo
      enddo

      !cont = number of faults to fix
      !tail = number of samples before the tail

      !identify the subfauts to freeze
      !set previous solution as prior
      !and change the weight (increase it)
      if ( green_mesh%rake_opt .eq. 1 ) then
       cont = 0
       k = 1
       do i=1,green_mesh%nsdip
        do j=1,green_mesh%nsstk
         if ( serie(k) .lt. timemark ) then
          !replace tail of time series for prior model
          ii=(k-1)*green_mesh%interp_i+1
          jj=k*green_mesh%interp_i
          !print *, 1, tail, ii, ii+tail-1
          green_mesh%p_model1d(1:tail,k) = green_mesh%model1(ii:ii+tail-1)
          weight(1:tail,k) = newweight
          cont = cont + 1
          idsub(cont) = k
         else
          !do not fix
         endif
         k = k + 1
        enddo
       enddo
      elseif ( green_mesh%rake_opt .eq. 2 ) then
       cont = 0
       k = 1
       do i=1,green_mesh%nsdip
        do j=1,green_mesh%nsstk
         if ( serie(k) .lt. timemark ) then
          !replace tail of time series for prior model
          ii=(k-1)*green_mesh%interp_i*2+1
          !print *, 1, 1, tail, k, ii, ii+tail-1
          green_mesh%p_model2d(1:tail,k) = green_mesh%model2(ii:ii+tail-1)
          !2nd component of the vector
          ii = ii + green_mesh%interp_i
          !print *, 2, 1, tail, k+green_mesh%msub, ii, ii+tail-1
          green_mesh%p_model2d(1:tail,k+green_mesh%msub) = green_mesh%model2(ii:ii+tail-1)
          weight(1:tail,k) = newweight
          cont = cont + 1
          idsub(cont) = k
         else
          !do not fix
         endif
         k = k + 1
        enddo
       enddo
      endif


      !Arrange weight back to the diagonal shape
      k=1
      do i=1,green_mesh%interp_i
       do j=1,green_mesh%msub
        green_mesh%diag(k) = weight(i,j)
        k=k+1
       enddo
      enddo


      write(*,*) ' Subfaults fixed:', cont

      open(63,file='weigth.dat',status='unknown')
      do i=1,green_mesh%msub*green_mesh%interp_i
       write(63,*) green_mesh%diag(i)
      enddo
      close(63)
      call write_prior(green_mesh)
      call write_model(green_mesh)

      !call prior_future(green_mesh,idsub,serie,cont,timewin)

      deallocate(serie,idsub,weight)
      endsubroutine identify_past



     subroutine prior_future(green_mesh,idsub,serie,cont,timewin)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      real,intent(inout):: timewin
      integer,intent(inout):: cont
      integer,intent(inout):: idsub(green_mesh%msub)
      real,intent(inout):: serie(green_mesh%msub)
      integer :: iunit, i, j, k, ii, jj, counter, counter2, counter3, counter4, counter5
      integer :: id_hypo, id_max, id_samp
      integer,dimension(:),allocatable :: board, samples
      real,dimension(:),allocatable :: valfar
      real :: farhypo, maxv

      allocate(board(green_mesh%msub))
      allocate(samples(green_mesh%msub))
      allocate(valfar(green_mesh%msub))
      valfar(:) = 0.
      board(:) = 0

      !need to identify the subfault nearest to the hypocenter
      id_hypo = 422

      !distance far from the hypocenter: set to half of maximum distance
      timewin = (real(green_mesh%samp4win(green_mesh%synwin-1))-1.)&
  &             *green_mesh%slipdt
      timewin = timewin*0.7
      ii=int(timewin/green_mesh%slipdt)+1
      timewin = (real(green_mesh%samp4win(green_mesh%synwin-1))-1.)&
  &             *green_mesh%slipdt
      jj=int(timewin/green_mesh%slipdt)+1
      print *, ' Sanmples as limit:', ii, jj


      counter=(id_hypo-1)*green_mesh%interp_i+1
      !read pre-estimated sample of the rupture
      iunit=22
      open(iunit,file=green_mesh%dat//'timespshift.dat',status='old',action='read')
      do i=1,green_mesh%msub
       read(iunit,*) samples(i)
       if ((( samples(i) .gt. ii ) .and. ( samples(i) .le. jj )) &
  &         .and. ( idsub(i) .ne. 0 )) then
        counter3=(i-1)*green_mesh%interp_i+(samples(i)-2)
        counter4=green_mesh%interp_i-(samples(i)-2)
        print *, samples(i), counter4
        green_mesh%model1(counter3:counter3+counter4) = &
  &          green_mesh%model1(counter:counter+counter4)
       endif
      enddo
      close(iunit)

     open(64,file='pred.mod1',status='unknown')
     do i=1,green_mesh%interp_i*green_mesh%msub
      write(64,*) green_mesh%model1(i)
     enddo
     close(64)


pause
      !Maximum slip rate from the close-to-bourder nodes
      !maxv = maxval(valfar)
      !do i=1,cont
      !  if ( maxv .eq. valfar(i) ) then
      !    id_max = idsub(i)
      !    counter = (idsub(i)-1)*green_mesh%interp_i+1
      !    do j=1,green_mesh%interp_i
      !     if ( maxv .eq. green_mesh%model1(counter) ) then
      !        id_samp = j
      !     endif
      !     counter = counter+1
      !    enddo
      !  endif
      !enddo

      id_max = 422
      !print *, id_max, 'id_max', id_samp

      !counter3 = (id_max-1)*green_mesh%interp_i+id_samp-2

      !distance far from the maximum slip-rate and out 
      !the active zone (future zone to break)
      !farhypo = (timewin*2.)*1000.
      do i=1,green_mesh%msub
        if ( (green_mesh%distcorr(i,id_max) .le. farhypo) .and. &
  &          (board(i) .ne. 1) ) then
         board(i) = 2
         !print *, board(i), i
         !counter = (i-1)*green_mesh%interp_i+samples(i)
         !counter2 = i*green_mesh%interp_i
         !counter4 = counter2-counter
         !green_mesh%model1(counter:counter2) != &
!  &              green_mesh%model1(counter3:counter3+counter4)
        endif
      enddo

     !open(64,file='pred.mod1',status='unknown')
     !do i=1,green_mesh%interp_i*green_mesh%msub
     ! write(64,*) green_mesh%model1(i)
     !enddo
     !close(64)

     deallocate(board,valfar,samples)
     endsubroutine prior_future


     subroutine detect_motion(green_mesh)

      !subroutine dedicated to the analysis of the slip-rate history
      !in order to detect the velocity and direction of the rupture
      !to furhter predict an initial guess of the future rupture

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer :: iunit, i, j, k, ii, jj, nsnap
      integer :: vel_x, vel_y, jump, ini, next_x, next_y, intersamp
      integer,dimension(:),allocatable :: snap, idmax_x, idmax_y, samples
      real,dimension(:),allocatable :: valfar
      integer,dimension(:,:),allocatable :: board, idcorr
      real,dimension(:,:),allocatable :: sliprate
      real :: farhypo, maxv, timeprev, interval
      integer :: node_bef, node_aft, substk, subdip, cont, cont2
      real :: corr_dist, idmax

      !interval to see the snapshots (seconds)
      interval = 1.
      intersamp = int(interval/green_mesh%slipdt) - 1

      !array to save snapshots
      allocate(valfar(green_mesh%msub))
      valfar(:)=0.

      !final time at this stage
      timeprev = real(green_mesh%samp4win(green_mesh%synwin-1))*green_mesh%slipdt

      allocate(samples(green_mesh%msub))
      iunit=22
      open(iunit,file=green_mesh%dat//'timespshift.dat',status='old',action='read')
      do i=1,green_mesh%msub
       read(iunit,*) samples(i)
      enddo


      !number of snapshots to look at      
      nsnap = (int(timeprev/green_mesh%slipdt)-2)/intersamp
      allocate(snap(nsnap),idmax_x(nsnap),idmax_y(nsnap))
      snap(:)=0
      idmax_x(:)=0
      idmax_y(:)=0

      ini = 2
      do ii=1,nsnap
       snap(ii) = ini + (ii-1)*intersamp

      k=1
      do i=1,green_mesh%nsdip
       do j=1,green_mesh%nsstk
        jump = (k-1)*green_mesh%interp_i+snap(ii)
        valfar(k) = green_mesh%model1(jump)
        k=k+1
       enddo
      enddo

     maxv = maxval(valfar)

      k=1
      do i=1,green_mesh%nsdip
       do j=1,green_mesh%nsstk
        if ( maxv .eq. valfar(k) ) then
         idmax_x(ii) = j
         idmax_y(ii) = i
        endif
        k=k+1
       enddo
      enddo
      print *, snap(ii), idmax_x(ii), idmax_y(ii)

     enddo

     vel_x = idmax_x(nsnap) - idmax_x(nsnap-1)
     vel_y = idmax_y(nsnap) - idmax_y(nsnap-1)
     vel_x = -2
     vel_y = -1
     next_x = idmax_x(nsnap) + vel_x
     next_y = idmax_y(nsnap) + vel_y
     node_bef = 422-(green_mesh%synwin-5)*36-(green_mesh%synwin-5)!(idmax_y(nsnap)-1)*green_mesh%nsstk+idmax_x(nsnap)
     node_aft = node_bef-(36+(green_mesh%synwin-5))!(next_y-1)*green_mesh%nsstk+next_x
     print *, node_bef, node_aft


     !ii = (node_bef-1)*green_mesh%interp_i+1
     !valfar(:) = 0.
     !valfar(1:green_mesh%interp_i) = green_mesh%model1(ii:ii+green_mesh%interp_i-1)
     !maxv = maxval(valfar)


     !Correlation distances
     corr_dist = 2.5*1000
     !Subfaults inside the correlation distances
     substk = int((corr_dist/1000.)*2./green_mesh%stk_s)
     subdip = int((corr_dist/1000.)*2./green_mesh%dip_s)
     !Array to save sliprates from the past
     allocate(sliprate(green_mesh%interp_i,substk*subdip))
     allocate(board(green_mesh%msub,2))   !past=1, future=2
     allocate(idcorr(substk*subdip,3))
     sliprate(:,:) = 0.
     board(:,:) = 0
     idcorr(:,:) = 0

     k=1
     cont=0
     cont2=0
     do i=1,green_mesh%nsdip
      do j=1,green_mesh%nsstk
       !Correlation in the past
       if ( green_mesh%distcorr(k,node_bef) .le. corr_dist) then
         cont=cont+1
         ii = (k-1)*green_mesh%interp_i+1
         idcorr(cont,1) = k
         board(k,1) = 1
print *, k, 'pas'
         sliprate(:,cont) = green_mesh%model1(ii:ii+green_mesh%interp_i-1)
       endif
       !Correlation in the future
       if ( green_mesh%distcorr(k,node_aft) .le. corr_dist) then
         cont2=cont2+1
         board(k,2) = 2
         idcorr(cont2,2) = k
print *, k, 'fut'
       endif
       k = k+1
      enddo
     enddo


     do i=1,cont2
      j = 1
      ii = 0
      do while ( ( ii .ne. 1) .and. (j .le. cont2) )
        if ( idcorr(i,2) .eq. idcorr(j,1)) then
          ii = 1
          !print *, idcorr(i,2), idcorr(j,1),i,j, 'no change'
        endif
        j=j+1
      enddo
      if ( ii .ne. 1 ) then
       jj = (idcorr(i,1)-1)*green_mesh%interp_i+1
       valfar(:) = 0.
       valfar(1:green_mesh%interp_i) = green_mesh%model1(jj:jj+green_mesh%interp_i-1)
       maxv = maxval(valfar)
       do k=1, green_mesh%interp_i
        if ( maxv .eq. valfar(k) ) then
         idmax = k
        endif
       enddo
       ii = samples(idcorr(i,2)) - 2
       j = green_mesh%interp_i - (ii-2)
       ii = (idcorr(i,2)-1)*green_mesh%interp_i+ii
       green_mesh%model1(ii:ii+j) = valfar(idmax-2:(idmax-2)+j)
      else 
       !print *, i, j-1, idcorr(i,2), idcorr(j-1,1), 'no change'
      endif
     enddo

     open(64,file='pred.mod1',status='unknown')
     do i=1,green_mesh%interp_i*green_mesh%msub
      write(64,*) green_mesh%model1(i)
     enddo
     close(64)


     !Smoothing the prediction
     !green_mesh%grad1(:) = green_mesh%model1(:)
     !call arrange_grad1d(green_mesh)
     !green_mesh%model1(:) = green_mesh%grad1(:)
     !green_mesh%grad1(:) = 0.

     open(65,file='pred.mod1s',status='unknown')
     do i=1,green_mesh%interp_i*green_mesh%msub
      write(65,*) green_mesh%model1(i)
     enddo
     close(65)


     deallocate(sliprate,board,idcorr,samples)
     deallocate(valfar,snap,idmax_x,idmax_y)
     endsubroutine detect_motion

