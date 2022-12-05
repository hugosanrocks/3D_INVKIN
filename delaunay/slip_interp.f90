      program slip_interp

       !COMMON VARIABLES
       IMPLICIT NONE
       integer ntime, ncomp, msub, stcomp, msubint
       real,dimension(:,:),allocatable :: tractionvec, green
       !=======================================!
       ! Forward model using Trapezoidal rule  !
       !=======================================!
       integer i, j, k, ii, jj, kk, cont, m
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
       stcomp = 120

       allocate(tractionvec(ntime,msub))


       !Size of subfaults along strike and dip
       write(*,*) 'dx coarse grid (km):' 
       read(*,*) dx
       dx = dx*1000.
       write(*,*) 'dy coarse grid (km):' 
       read(*,*) dy
       dy = dy*1000.
       write(*,*) 'nx and ny coarse grid (subfaults):' 
       read(*,*) nx, ny

       call read_slip(tractionvec,msub,ntime)


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

       write(*,*) 'Read grid for interpolation:'
       read(*,*) msubint
       allocate(vxint(msubint),vyint(msubint))
       allocate(green(ntime,msubint))
       !Read grid
       open(iunit,file='grid.in',status='old',action='read')
       do i=1,msubint
        read(iunit,*) vxint(i), vyint(i)
       enddo
       close(iunit)
       print *, vxint(1:4)
       print *, vyint(1:4)

       !scalfac=green_mesh%moment/(green_mesh%msu
       allocate(fault(ntime,nx,ny))
      ! allocate(green(green_mesh%trac_i,green_mesh%stcomp*green_mesh%ncomp*nxint*nyint))
       allocate(fault1(nx,ny))

       do i=1,ntime
         k=1
         do ii=1,ny
          do jj=1,nx
           fault(i,jj,ii) = tractionvec(i,k)
           k = k + 1
          enddo
         enddo
         fault1(:,:) = fault(i,:,:)
         k = 1
          do jj=1,msubint
            xwant=vxint(jj)
            ywant=vyint(jj)
write(*,*) xwant, ywant
            call FIND2D(xwant,ywant,res,vx,vy,fault1,nx,ny)
            green(i,k) = res
            k = k + 1
           enddo
       enddo

       !Divide green functions into 2 records to write them
       iunit=16
       OPEN(iunit,file='model_interp.out',status='unknown',&
 &          form='formatted',action='write')
         do i=1,ntime
          write(iunit,*) green(i,:)
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
 
     endprogram slip_interp




      subroutine read_slip(tractionvec,msub,ntime)

       !COMMON VARIABLES
       IMPLICIT NONE
       integer iunit, i
       integer,intent(inout) :: msub, ntime
       real,intent(inout) :: tractionvec(ntime,msub)

      
       iunit=16
        write(*,*) ' Reading slip-rate functions from dat/vitesse.out '
            OPEN(iunit,file='model_target.dat',status='old',&
 &          form='formatted',action='read')
           do i=1,ntime
            read(iunit,*) tractionvec(i,:)
           enddo
       close(iunit)


      end subroutine read_slip



