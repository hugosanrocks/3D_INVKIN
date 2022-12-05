      program green_interp

       !COMMON VARIABLES
       IMPLICIT NONE
       integer ntime, nsta, ncomp, msub, stcomp, msubint
       real,dimension(:,:),allocatable :: tractionvec, green
       !=======================================!
       ! Forward model using Trapezoidal rule  !
       !=======================================!
       integer i, j, k, ii, jj, kk, cont, m, im, jm
       integer fi, fj, fk, count2, count1, iunit, reclent
       real scalfac
       real, dimension(:,:,:), allocatable :: fault
       real, dimension(:,:), allocatable :: fault1

       !Variables for 2D integration
       real,dimension(:),allocatable :: vx, vy, vxint, vyint
       real dx, dy, res, xwant, ywant
       integer nx, ny, nxint, nyint, col1, col2, col3, col4

       write(*,*) 'Time samples in Green functions:'
       read(*,*) ntime
       write(*,*) 'Subfaults msub:'
       read(*,*) msub
       ncomp = 3
       nsta = 40
       stcomp = nsta*ncomp

       allocate(tractionvec(ntime,msub*ncomp*stcomp))


       !Size of subfaults along strike and dip
       write(*,*) 'dx coarse grid (km):' 
       read(*,*) dx
       dx = dx*1000.
       write(*,*) 'dy coarse grid (km):' 
       read(*,*) dy
       dy = dy*1000.
       write(*,*) 'nx and ny coarse grid (subfaults):' 
       read(*,*) nx, ny

       call read_time(tractionvec,nx,ny,msub,ncomp,stcomp,ntime)


       allocate(vx(nx),vy(ny))

       !Build coordinate vectors for 2D integration
       vx(1) = 0.*1000.
       vy(1) = 0.*1000.
       do i=2,nx
         vx(i) = vx(i-1)+dx
       enddo
       do i=2,ny
         vy(i) = vy(i-1)+dy
       enddo
       print *, vx(1:4)
       print *, vy(1:4)
       print *, vx(nx-4:nx)
       print *, vy(ny-4:ny)

       write(*,*) 'reading random grid (subfaults):'
       read(*,*) msubint
       allocate(vxint(msubint),vyint(msubint))
       allocate(green(ntime,stcomp*ncomp*msubint))
       !Build X and Y vectors for 2D interpolation
       write(*,*) 'reading grid form file'
       iunit=11
       open(iunit,file='grid.in',status='old',action='read')
       do i=1,msubint
        read(iunit,*) vxint(i), vyint(i)
       enddo
       vxint(:) = vxint(:)*1000.
       vyint(:) = vyint(:)*1000.
       close(iunit)
       print *, vxint(1:4)
       print *, vyint(1:4)

       allocate(fault(ntime,nx,ny))
       allocate(fault1(nx,ny))

      do i=1,nsta
      print *, '***sta', i
       do j=1,ncomp
        do m=1,ncomp
         k=1
         do ii=1,ny
          do jj=1,nx
           kk = m + (i-1)*ncomp**2 & 
  &             + (k-1)*stcomp*ncomp &
  &             + (j-1)*ncomp
           !used to check elements
           !print *, 'conv t', kk, 'f',k
           fault(:,jj,ii) = tractionvec(:,kk)
           k = k + 1
          enddo
         enddo
         do cont = 1,ntime
          fault1(:,:) = fault(cont,:,:)
          if (( cont .eq. 10 ) .and. ( i .eq. 10 )) then
            do jm=1,ny
             do im=1,nx
               write(80,*) fault1(im,jm)
             enddo
            enddo
          endif
          k = 1
          do ii=1,msubint
            kk = m + (i-1)*ncomp**2 &
  &             + (k-1)*stcomp*ncomp &
  &             + (j-1)*ncomp
            xwant=vxint(ii)
            ywant=vyint(ii)
            call FIND2D(xwant,ywant,res,vx,vy,fault1,nx,ny)
            green(cont,kk) = res
            k = k + 1
          enddo
          if (( cont .eq. 10 ) .and. ( i .eq. 10 )) then
             do im=1,msubint
               write(81,*) green(cont,im)
             enddo
             !stop
          endif
         enddo
        enddo
       enddo
      enddo

       !Divide green functions into 12 records to write them
       col1=msubint*stcomp*ncomp/12
       reclent=ntime*col1
       iunit=16
       OPEN(iunit,file='TRACT_time_rand.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
         do i=1,12
          write(iunit,rec=i) green(:,1+((i-1)*col1):col1*i)
         enddo
       close(iunit)

       !To check, read the green functions again
       !green(:,:)=0.

       !iunit=16
       !OPEN(iunit,file='dat/TRACT_time_interp.bin',status='unknown',&
! &          form='unformatted',ACCESS='DIRECT',recl=reclent)
       !  read(iunit,rec=1) green(:,1:col1)
       !  read(iunit,rec=2) green(:,col1+1:col2)
       !close(iunit)

       !do i=1,ntime
       ! write(77,*) green(i,1), green(i,col1+1)
       !enddo


      deallocate(fault,fault1)
      deallocate(vx,vy,vxint,vyint)
      deallocate(green)
 
     endprogram green_interp




      subroutine read_time(tractionvec,nsstk,nsdip,msub,ncomp,stcomp,ntime)

       !COMMON VARIABLES
       IMPLICIT NONE
       integer bit, iunit, i, j, k, l, reclent
       integer,intent(inout) :: nsstk, nsdip, msub, stcomp, ncomp, ntime
       real,intent(inout) :: tractionvec(ntime,msub*stcomp*ncomp)
       integer col1
       integer ini, fin, divide

       divide = 12   !packet size divided

       bit=1
       ! Total elements used for traction matrix (frequency)
       col1=nsstk*nsdip*stcomp*ncomp/divide
       reclent=ntime*col1

       iunit=16
        write(*,*) ' Reading Green functions from dat/TRACT_time.bin '
            OPEN(iunit,file='TRACT_time.bin',status='old',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
        do i=1,divide
            ini = 1+(col1*(i-1))
            fin = col1*i
            write(*,*) ' Read Green s functions: ',i,'/',12
            read(iunit,rec=i) tractionvec(:,ini:fin)
        enddo
       close(iunit)




      end subroutine read_time




