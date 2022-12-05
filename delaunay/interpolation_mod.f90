
      module interpolate_grid

      contains

      subroutine interpolate_green_ref(green_mesh,flag_read)

      !========================================================!
      !          Subroutine interpolate_green_ref              !
      !                                                        !
      ! Reads info from the Green functions of the reference   !
      ! grid and then it interpolates it into the unstructured !
      ! one.                                                   !
      !========================================================!

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh
       integer, intent(inout) :: flag_read

       write(*,*) '================================================='
       write(*,*) ' An unstructured grid (dat/nodes.in) will be used'
       write(*,*) '================================================='

       ! Read the info from reference grid
       if ( flag_read .eq. 1 ) then
         ! Read info about reference grid
         call read_green_ref_info(green_mesh)
         ! Interpolate reference grid to the unstructured grid
         call green_interpol(green_mesh%trac_i,green_mesh%nsta,&
  &                          green_mesh%ncomp,green_mesh%msub,&
  &                          green_mesh%nnodes,green_mesh%nsstk,&
  &                          green_mesh%nsdip,green_mesh%stk_ref,&
  &                          green_mesh%dip_ref,green_mesh%nodes_coor,&
  &                          green_mesh%tractionvec_ref,green_mesh%tractionvec)
       else
         ! Interpolate reference grid to the unstructured grid
         call green_interpol(green_mesh%trac_i,green_mesh%nsta,&
  &                          green_mesh%ncomp,green_mesh%msub,&
  &                          green_mesh%nnodes,green_mesh%nsstk,&
  &                          green_mesh%nsdip,green_mesh%stk_ref,&
  &                          green_mesh%dip_ref,green_mesh%nodes_coor,&
  &                          green_mesh%tractionvec_ref,green_mesh%tractionvec)
       endif

       ! Back to correct number of nodes
       green_mesh%msub = green_mesh%nnodes

       write(*,*) '================================================='
       write(*,*) ' Ready Green functions for unstructured grid     '
       write(*,*) '================================================='

      end subroutine interpolate_green_ref


      subroutine read_green_ref_info(green_mesh)

      !========================================================!
      !          Subroutine read_green_info                    !
      !                                                        !
      ! Reads info of the cartesian grid that has to be        !
      ! interpolated to the unstructured grid.                 !
      !                                                        !
      !========================================================!

       !COMMON VARIABLES
       use delaunay
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh
       integer i, iunit

       !Read info of underlying regular grid
       write(*,*) ' Reading info of underlying cartesian grid '
       iunit = 22
       open(iunit,file='dat/green_interp.info',status='old',&
  &         action='read')
       read(iunit,*) green_mesh%trac_i, green_mesh%msub
       write(*,*) ' Time samples in Green functions: ', green_mesh%trac_i
       write(*,*) ' Number of nodes in regular grid: ', green_mesh%msub
       read(iunit,*) green_mesh%ncomp, green_mesh%nsta
       write(*,*) ' Number of stations: ', green_mesh%nsta
       write(*,*) ' Components of seismograms: ', green_mesh%ncomp

       !Size of subfaults along strike and dip
       read(iunit,*) green_mesh%dstk, green_mesh%ddip
       write(*,*) ' dx, dy underlying grid (m): ', green_mesh%dstk, green_mesh%ddip
       read(iunit,*) green_mesh%nsstk, green_mesh%nsdip
       write(*,*) ' nx, ny underlying grid (nodes): ', green_mesh%nsstk, green_mesh%nsdip

       call read_green_ref(green_mesh%tractionvec_ref,green_mesh%nsstk,&
  &                        green_mesh%nsdip,green_mesh%nsta,green_mesh%ncomp,&
  &                        green_mesh%trac_i)

       !Build coordinate vectors for 2D integration
       green_mesh%stk_ref(1) = 0.
       green_mesh%dip_ref(1) = 0.
       do i=2,green_mesh%nsstk
         green_mesh%stk_ref(i) = green_mesh%stk_ref(i-1)+green_mesh%dstk
       enddo
       do i=2,green_mesh%nsdip
         green_mesh%dip_ref(i) = green_mesh%dip_ref(i-1)+green_mesh%ddip
       enddo
       write(*,*) '================================================='
       write(*,*) ' Grid of reference:                              '
       write(*,*) '================================================='
       write(*,*) ' vx(1:4):  ', green_mesh%stk_ref(1:4)
       write(*,*) ' vy(1:4):  ', green_mesh%dip_ref(1:4)
       write(*,*) ' vx(nx-3): ', green_mesh%stk_ref(green_mesh%nsstk-3:green_mesh%nsstk)
       write(*,*) ' vy(ny-3): ', green_mesh%dip_ref(green_mesh%nsdip-3:green_mesh%nsdip)
       write(*,*) ' Reading random grid (dat/nodes.in) '
       read(iunit,*) green_mesh%nnodes
       close(iunit)

       !Build X and Y vectors for 2D interpolation
       call read_nodes(green_mesh%nnodes,green_mesh%nodes_coor)

       write(*,*) '================================================='
       write(*,*) ' Unstructured grid:                              '
       write(*,*) '================================================='
       write(*,*) ' vxint(1:4):  ', green_mesh%nodes_coor(1,1:4)
       write(*,*) ' vyint(1:4):  ', green_mesh%nodes_coor(2,1:4)
       write(*,*) ' vxint(msubint-3): ', green_mesh%nodes_coor(1,green_mesh%nnodes-3:green_mesh%nnodes)
       write(*,*) ' vyint(msubint-3): ', green_mesh%nodes_coor(2,green_mesh%nnodes-3:green_mesh%nnodes)


      endsubroutine read_green_ref_info



      !========================================================!
      !          Subroutine green_interpol                     !
      !                                                        !
      ! Performs the Green function interpolation according    !
      ! to all the given parameters of the cartesian and non-  !
      ! structure grids.                                       !
      !                                                        !
      !========================================================!

      subroutine green_interpol(ntime,nsta,ncomp,msub,msubint,&
  &                             nx,ny,vx,vy,pos_int,&
  &                             tractionvec,green)

       !COMMON VARIABLES
       use interplib
       IMPLICIT NONE
       integer, intent(inout) :: ntime, nsta, ncomp, msub, msubint
       integer, intent(inout) :: nx, ny
       real, intent(inout) :: tractionvec(ntime,msub*ncomp*nsta*ncomp)
       real, intent(inout) :: green(ntime,nsta*ncomp*ncomp*msubint)
       real, intent(inout) :: vx(nx), vy(ny), pos_int(2,msubint)
       !Counters
       integer i, j, k, ii, jj, kk, cont, m
       !For interpolation
       real, dimension(:,:,:), allocatable :: fault
       real, dimension(:,:), allocatable :: fault1
       real, dimension(:), allocatable :: vxint, vyint
       real res, xwant, ywant

       allocate(fault(ntime,nx,ny))
       allocate(fault1(nx,ny))
       allocate(vxint(msubint),vyint(msubint))

       vxint(:) = pos_int(1,:)
       vyint(:) = pos_int(2,:)

       do i=1,nsta
        do j=1,ncomp
         do m=1,ncomp
          k=1
          do ii=1,ny
           do jj=1,nx
            kk = m + (i-1)*ncomp**2 & 
  &             + (k-1)*nsta*ncomp*ncomp &
  &             + (j-1)*ncomp
            !used to check elements
            !print *, 'conv t', kk, 'f',k
            fault(:,jj,ii) = tractionvec(:,kk)
            k = k + 1
           enddo
          enddo
          do cont = 1,ntime
           fault1(:,:) = fault(cont,:,:)
           k = 1
           do ii=1,msubint
             kk = m + (i-1)*ncomp*ncomp &
  &             + (k-1)*nsta*ncomp*ncomp &
  &             + (j-1)*ncomp
             xwant=vxint(ii)
             ywant=vyint(ii)
             call FIND2D(xwant,ywant,res,vx,vy,fault1,nx,ny)
             green(cont,kk) = res
             k = k + 1
           enddo
          enddo
         enddo
        enddo
       enddo

       ! Write the file with Green function on the unstructured grid
       call write_green_rand(green,msubint,nsta,ncomp,ntime)


      !Flush memory
      deallocate(fault,fault1,vxint,vyint)
      endsubroutine green_interpol


      !========================================================!
      !          Subroutine read_green                         !
      !                                                        !
      ! Reads the Green functions from the cartesian grid      !
      ! that have to be interpolated for the non-structured    !
      ! grid.                                                  !
      !                                                        !
      !========================================================!

      subroutine read_green_ref(tractionvec,nsstk,nsdip,nsta,ncomp,ntime)

       !COMMON VARIABLES
       IMPLICIT NONE
       integer bit, iunit, i, reclent
       integer,intent(inout) :: nsstk, nsdip, nsta, ncomp, ntime
       real,intent(inout) :: tractionvec(ntime,nsstk*nsdip*nsta*ncomp*ncomp)
       integer col1
       integer ini, fin, divide

       divide = 12   !packet size divided

       bit=4
       ! Total elements used for traction matrix (frequency)
       col1=nsstk*nsdip*nsta*ncomp*ncomp/divide
       reclent=bit*ntime*col1

       iunit=16
        write(6, *) ' Reading Green functions from dat/TRACT_time.bin '
            OPEN(iunit,file='dat/TRACT_time.bin',status='old',&
 &          form='unformatted',ACCESS='DIRECT',action='read',&
 &          recl=reclent)
        do i=1,divide
            ini = 1+(col1*(i-1))
            fin = col1*i
            write(6, *) ' Read Green s functions: ',i,'/', divide
            read(iunit,rec=i) tractionvec(:,ini:fin)
        enddo
       close(iunit)

      end subroutine read_green_ref


      !========================================================!
      !          Subroutine write_green                         !
      !                                                        !
      ! Reads the Green functions from the cartesian grid      !
      ! that have to be interpolated for the non-structured    !
      ! grid.                                                  !
      !                                                        !
      !========================================================!

      subroutine write_green_rand(tractionvec,nnodes,nsta,ncomp,ntime)

       !COMMON VARIABLES
       IMPLICIT NONE
       integer bit, iunit, i, reclent
       integer,intent(inout) :: nnodes, nsta, ncomp, ntime
       real,intent(inout) :: tractionvec(ntime,nnodes*nsta*ncomp*ncomp)
       integer col1
       integer ini, fin, divide

       divide = 12   !packet size divided

       bit=4
       ! Total elements used for traction matrix (frequency)
       col1=nnodes*nsta*ncomp*ncomp/divide
       reclent=bit*ntime*col1

       iunit=16
        write(*,*) ' Writing Green functions at dat/TRACT_time_rand.bin '
            OPEN(iunit,file='dat/TRACT_time_rand.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',action='write',&
 &          recl=reclent)
        do i=1,divide
            ini = 1+(col1*(i-1))
            fin = col1*i
            write(6, *) ' Writing Green s functions: ',i,'/', divide
            write(iunit,rec=i) tractionvec(:,ini:fin)
        enddo
       close(iunit)

      end subroutine write_green_rand


      !========================================================!
      !          Subroutine sliprate_interpol()                    !
      !                                                        !
      ! Performs the sliprate function interpolation according !
      ! to all the given parameters of the cartesian and non-  !
      ! structure grids.                                       !
      !                                                        !
      !========================================================!
      subroutine sliprate_interpol(nodesin,gridin,nelem,corners,nodesout,gridout,ntime,slipratein,sliprateout)

       !COMMON VARIABLES
       use triangles
       IMPLICIT NONE
       integer ntime, nodesin, nodesout, nelem
       real :: slipratein(ntime,nodesin), sliprateout(ntime,nodesout)
       real :: gridin(nodesin,2), gridout(nodesout,2)
       integer :: corners(nelem,3)
       !Variables for 2D integration
       real :: triang(3,2)
       real :: val(3), val_int
       integer, dimension(:), allocatable :: id
       !For interpolation
       integer :: i, j, k, check
       real :: xwant, ywant

       allocate(id(nodesout))
       !call reshape_sliprate(slipratein,ntime,nodesin)

       do i=1,nodesout
        xwant = gridout(i,1)
        ywant = gridout(i,2)
        check = 0
        j = 1
        do while ( (check .eq. 0) .and. (j .le. nelem) )
         triang(1,:) = [gridin(corners(j,1),1), gridin(corners(j,1),2)]
         triang(2,:) = [gridin(corners(j,2),1), gridin(corners(j,2),2)]
         triang(3,:) = [gridin(corners(j,3),1), gridin(corners(j,3),2)]
         check = point_in_triang(xwant,ywant,triang(1,1),triang(1,2),&
  &              triang(2,1),triang(2,2),triang(3,1),triang(3,2))
         j = j + 1
        enddo
        id(i) = j-1
        do k=1,ntime
         val(1:3) = [slipratein(k,corners(j-1,1)), &
  &                 slipratein(k,corners(j-1,2)), &
  &                 slipratein(k,corners(j-1,3))]
         val_int = triangle_interpol(xwant,ywant,triang(:,1),triang(:,2),val)
         sliprateout(k,i) = val_int
        enddo
       enddo
 
      deallocate(id)
      endsubroutine sliprate_interpol


      !subroutine reshape_sliprate(tractionvec,msub,ntime)

      ! !COMMON VARIABLES
      ! IMPLICIT NONE
      ! integer iunit, i
      ! integer,intent(inout) :: msub, ntime
      ! real,intent(inout) :: tractionvec(ntime,msub)
      

      !endsubroutine reshape_sliprate

      endmodule interpolate_grid
