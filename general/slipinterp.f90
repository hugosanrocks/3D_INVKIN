       program slipinterp


       implicit none
       integer siz, i, j, k, iunit, nsubf, nt
       integer sizint, nsubfint
       integer inx, iny, nxf, nyf, nxint, nyint
       real t, dt
       real, dimension(:,:), allocatable :: coor, coorint, slip, slipint
       real, dimension(:,:), allocatable :: snap, snapint
       real, dimension(:), allocatable :: xf, yf, xint, yint
       real xwant, ywant, value
       real dx, dy, x, y

       nsubf = 288
       nxf = 24
       nyf = 12
       nt = 225
       dt = 0.04
       nsubfint = 648
       nxint = 36
       nyint = 18

       iunit = 6

       !Read the slip solution
       open(iunit,file='modelpri.dat',status='old',action='read')
       siz = nsubf * nt
       sizint = nsubfint * nt
       allocate(slip(3,siz),slipint(3,sizint))
       read(iunit,*) slip
       close(iunit)

       allocate(coor(3,nsubf),coorint(3,nsubfint))

       !Read the coordinates of solution
       open(iunit,file='dat/fault.dat',status='old',action='read')
       read(iunit,*) coor
       close(iunit)

       !Read coordinates where to interpolate
       open(iunit,file='dat/fault.dat2',status='old',action='read')
       read(iunit,*) coorint
       close(iunit)



       allocate(snap(nxf,nyf),snapint(nxint,nyint))
       allocate(xf(nxf),yf(nyf),xint(nxint),yint(nyint))

       !Build X and Y vectors for 2D interpolation
       dx = 1.5
       dy = 1.5
       x = 0.5
       do i = 1, nxf
         xf(i) = x
         x = x + dx
       enddo
       y = 0.5
       do i = 1, nyf
         yf(i) = y
         y = y + dy
       enddo
       print *, xf(nxf)

       !Build X and Y vectors for 2D interpolation
       dx = 1.
       dy = 1.
       x = 0.5
       do i = 1, nxint
         xint(i) = x
         x = x + dx
       enddo
       y = 0.5
       do i = 1, nyint
         yint(i) = y
         y = y + dy
       enddo
       print *, xint(nxint)



       !Interpolation at each snapshot
       do j = 1 , 3                     !each component
        do i = 1 , nt
        !Arrange slip as snapshots
        k = i
         do iny = 1, nyf
          do inx = 1,nxf               !each subfault
            snap(inx,iny) = slip(j,k)
            k = k + nt
          enddo
         enddo
        !===================================================!
        !Interpolate, arrange as snapshot and save as vector
        k = i
         do iny = 1, nyint
          do inx = 1, nxint
           xwant = xint(inx)
           ywant = yint(iny)
           call FIND2D(xwant,ywant,value,xf,yf,snap,nxf,nyf)
           snapint(inx,iny) = value
           slipint(j,k) = snapint(inx,iny)
           k = k + nt
          enddo
         enddo
!         write(44,*) snapint(28,12)
         if (i .eq. 40) then
         do iny = 1, nyint
          do inx = 1,nxint               !each subfault
          !  write(44,*) snapint(inx,iny)
          enddo
         enddo
         endif
         !---------------------------------------------------!
        enddo
       enddo


       !Read the slip solution
       open(iunit,file='dat/modelinterp.dat',status='unknown')
       do i = 1, sizint
       write(iunit,*) slipint(:,i)
       enddo
       close(iunit)



       deallocate(xf,yf,xint,yint)
       deallocate(snap,snapint)
       deallocate(slip,slipint)
       deallocate(coor,coorint)
       endprogram slipinterp
