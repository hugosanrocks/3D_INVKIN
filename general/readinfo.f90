      !==========================================================!
      !        Subroutine read_info                              !
      !                                                          !
      ! This subroutine is in charge to read all the information !
      ! necessary to run the kinematic source inversion          !
      ! The following files are read by this subroutine          !
      !
      ! dat/simul.info:     Information about stress-state tensors!
      ! dat/focal.info:     Local focal mechanism information     !
      ! dat/sliprate.info:  Source time history informaiton       !
      ! dat/fwioption.info: Choice of optimization strategy       !
      ! dat/syn.info:       Information about forward modeling    !
      ! dat/filter.info:    Choice of time filtering              !
      ! dat/green.info:     Information about underlying mesh     !
      !===========================================================!

      subroutine read_info(green_mesh)

      IMPLICIT NONE
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer iunit, i
      REAL pi/3.1415926535897932/

      !only used for I/O issues
      REAL :: vector2(2), vector3(3)

      !===============================================================!
      ! Setting input and output directories                          !
      !===============================================================!
      ! Directory containing input data files
      green_mesh%dat='dat/'
      ! Directory containing output data files
      green_mesh%out='out/'
      !===============================================================!

      !Create files for saving data and model norms
      iunit=44
      open(iunit,file=green_mesh%out//'data_norm.out',status='unknown')
      close(iunit)
      open(iunit,file=green_mesh%out//'model_norm.out',status='unknown')
      close(iunit)

      !create file to monitor allocations
      open(iunit,file=' alloc_monitor.out',status='unknown',action='write')
       write(iunit,*) ' ########### Allocation monitor ############ '
      close(iunit)



      WRITE(6, *) '================================================='
      write(6,*) ' Reading input information '
      WRITE(6, *) '================================================='

      !===================================================================!
      !           Reading simul.info                                      !
      !===================================================================!
      write(6,*) ' Reading Green s functions '
      open(iunit,file=green_mesh%dat//'simul.info',&
  &        status='old',action='read')
       read(iunit,*) green_mesh%rt, green_mesh%ot
       read(iunit,*) green_mesh%nsta, green_mesh%ncomp, green_mesh%msub
       read(iunit,*) green_mesh%simsam, green_mesh%simt, green_mesh%simdt
       read(iunit,*) green_mesh%stress_opt
       read(iunit,*) green_mesh%interp_i, green_mesh%slipdt, &
  &                  green_mesh%ntint, green_mesh%trac_i, &
  &                  green_mesh%lenobs 
      close(iunit)
      !===================================================================!


      !===================================================================!
      !           Reading focal.info                                      !
      !===================================================================!

      !======Rake geometry=================!

      !Relative motion of hanging wall
             !^ 270ª dip -
             !=================!> 0ª strike -
     !180ª+< !                 !
     !strike !                 !
             !=================!

      write(6,*) ' Reading focal mechanism '
      open(iunit,file=green_mesh%dat//'focal.info',&
  &        status='old',action='read')
      read(iunit,*) green_mesh%dir_n
      allocate(green_mesh%rake_i(green_mesh%dir_n))
      call monitor_alloc('54%rake_i ',' allo')
      allocate(green_mesh%vnorm(green_mesh%dir_n,3))
      call monitor_alloc('55%vnorm  ',' allo')
      allocate(green_mesh%vstk(green_mesh%dir_n,3))
      call monitor_alloc('56%vstk   ',' allo')
      allocate(green_mesh%vdip(green_mesh%dir_n,3))
      call monitor_alloc('57%vdip   ',' allo')
      allocate(green_mesh%vslip(green_mesh%dir_n,3))
      call monitor_alloc('58%vslip  ',' allo')
      allocate(green_mesh%vslip2(green_mesh%dir_n,2))
      call monitor_alloc('59%vslip2 ',' allo')
      allocate(green_mesh%slipm(green_mesh%dir_n,2,3))
      call monitor_alloc('60%slipm  ',' allo')
      allocate(green_mesh%vnorm_i(green_mesh%msub))
      call monitor_alloc('61%vnorm_i',' allo')
      allocate(green_mesh%rake_lim(green_mesh%dir_n,3))
      call monitor_alloc('62%rake_li',' allo')
      allocate(green_mesh%stk_i(green_mesh%dir_n))
      call monitor_alloc('63%stk_i  ',' allo')
      allocate(green_mesh%dip_i(green_mesh%dir_n))
      call monitor_alloc('64%dip_i  ',' allo')
      do i = 1, green_mesh%dir_n
       green_mesh%dir_i = i
       read(iunit,*) green_mesh%stk, green_mesh%dip, green_mesh%rak
       green_mesh%stk_i(i) = green_mesh%stk
       green_mesh%dip_i(i) = green_mesh%dip
       green_mesh%rake_i(i) = green_mesh%rak
       !Convert angles to radiants
       write(6,*) ' Focal mechanism', i, green_mesh%stk, green_mesh%dip, green_mesh%rak
       green_mesh%stk=green_mesh%stk*(2.0*pi/360.0)
       green_mesh%dip=green_mesh%dip*(2.0*pi/360.0)
       green_mesh%rak=(green_mesh%rak*(2.0*pi/360.0))
       call vectors_geometry(green_mesh)
      enddo
      !read(iunit,*) green_mesh%stk, green_mesh%dip, green_mesh%rak
      read(iunit,*) green_mesh%lenstk, green_mesh%lendip
      read(iunit,*) green_mesh%nsstk, green_mesh%nsdip
      read(iunit,*) green_mesh%stk_s, green_mesh%dip_s
      read(iunit,*) green_mesh%rake_opt
      read(iunit,*) green_mesh%hypo(:)
      !read(iunit,*) green_mesh%vrup
      read(iunit,*) green_mesh%rak_case
      do i = 1,green_mesh%dir_n
       read(iunit,*) vector2(1:2)
       green_mesh%rake_lim(i,1:2) = vector2(1:2)
       write(6, *) ' Rake limits:', vector2(1:2)
      enddo
      close(iunit)
      green_mesh%msub=green_mesh%nsstk*green_mesh%nsdip

      !===================================================================!


      !===================================================================!
      !           Reading sliprate.info                                      !
      !===================================================================!

      ! Read slip rate point source information
      write(6,*) ' Reading inital model '
      open(iunit,file=green_mesh%dat//'sliprate.info',status='old',action='read')
       read(iunit,*) green_mesh%slipsam, green_mesh%slipt
       read(iunit,*) green_mesh%moment
      close(iunit)

      !dt for slip-rate interpolation in time
      green_mesh%intdt = green_mesh%slipt &
  &                      / (real(green_mesh%slipsam)-1.)
      write(6, *) ' dt of slip-rates not interpolated: ', green_mesh%slipdt
      write(6, *) ' dt for slip-rate interpolation: ', green_mesh%intdt

      IF (green_mesh%debug .eqv. .true.) THEN
        ! Write simulation information to check
        green_mesh%debug_i = 1
        call write_debug(green_mesh)
      ENDIF

      !===================================================================!



      !===================================================================!
      !           Reading fwioption.info                                      !
      !===================================================================!

      !Read the option selected
      write(6,*) ' Reading optimization options '
      open(iunit,file=green_mesh%dat//'fwioption.info',status='unknown')
      read(iunit,*) green_mesh%fwiopt
      read(iunit,*) green_mesh%niter_max
      read(iunit,*) green_mesh%weig
      read(iunit,*) green_mesh%lam1, green_mesh%lam2, green_mesh%lam3, green_mesh%lam4, green_mesh%lam5
      read(iunit,*) green_mesh%quota1, green_mesh%quota2, green_mesh%quota3, green_mesh%quota4
      read(iunit,*) green_mesh%lb, green_mesh%ub
      read(iunit,*) green_mesh%depth_opt, green_mesh%depth_coef
      read(iunit,*) green_mesh%percent_win
      read(iunit,*) green_mesh%percent_pri
      read(iunit,*) green_mesh%new_weight
      read(iunit,*) green_mesh%wintochange
      close(iunit)
      !read map of focal mechanism
      call detect_focmec(green_mesh)
      !===================================================================!

      !===================================================================!
      !           Reading syn.info                                      !
      !===================================================================!
      !Write information about synthetics
      write(6,*) ' Reading inforamtion for synthetics '
      OPEN(iunit,file=green_mesh%dat//'syn.info',status='unknown')
      read(iunit,*) green_mesh%for_opt                    !null or prior information
      read(iunit,*) green_mesh%optm, green_mesh%mext, green_mesh%mxint     !extend model length
      read(iunit,*) green_mesh%wininv, green_mesh%dowin  !number of time window to forward
      read(iunit,*) green_mesh%synwin
      allocate(green_mesh%samp4win(green_mesh%wininv))
      call monitor_alloc('65%samp4wi',' allo')
      read(iunit,*) green_mesh%samp4win(:)
      close(iunit)
      !===================================================================!


      !===================================================================!
      !           Reading filter.info                                      !
      !===================================================================!
      open(iunit,file=green_mesh%dat//'filter.info',&
  &        status='old',action='read')
      read(iunit,*) green_mesh%f_obs, green_mesh%f_syn, green_mesh%f_grad
      read(iunit,*) green_mesh%fc_obs, green_mesh%fc_syn, green_mesh%fc_grad
      read(iunit,*) green_mesh%order
      read(iunit,*) green_mesh%filtgrad
      read(iunit,*) green_mesh%r1, green_mesh%r2, green_mesh%repeat
      read(iunit,*) green_mesh%lap_opt
      read(iunit,*) green_mesh%xhyp, green_mesh%zhyp
      close(iunit)
      !===================================================================!

      !===================================================================!
      !           Reading green.info                                      !
      !===================================================================!
      !input related to mesh
       open(iunit,file=green_mesh%dat//'green.info',status='old',&
  &         action='read')
       !optional mesh, regular or iregular
       read(iunit,*) green_mesh%opt_mesh
       !nodes along strike and dip of underlying mesh
       read(iunit,*) green_mesh%nsstk, green_mesh%nsdip
       !if iregular grid used, choose gaussian points
        if ( green_mesh%opt_mesh .eq. 3 ) then
         read(iunit,*) green_mesh%n_gauss, green_mesh%siz_gauss
         read(iunit,*) green_mesh%nnodes
         write(6,*) ' Triangular mesh is on '
         write(6,*) ' # of nodes in the grid: ', green_mesh%nnodes
        endif
      close(iunit)
      !===================================================================!


      WRITE(6, *) '================================================='
      write(6,*) ' Finish reading input information '
      WRITE(6, *) '================================================='

      end subroutine read_info




      !============================================================!
      ! Subroutine:     triangle_elements
      ! Given a mesh of nodes, a delaunay triangulation is
      ! performed. The id of the nodes forming the triangular
      ! elements are saved in  while all the nodes are saved in
      ! 
      !============================================================!

      subroutine read_triang_elem(green_mesh)

      use delaunay
      IMPLICIT NONE
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer triangle_order
      integer, dimension(:,:), allocatable :: triangle_neighbor

      ! Read node coordinates from dat/nodes.in
      call read_nodes(green_mesh%nnodes,green_mesh%nodes_coor)

      ! Determine the Delaunay triangulation.
      triangle_order = 3          !three nodes per triangle

      ! Allocate memory
      allocate(triangle_neighbor(triangle_order,triangle_order*green_mesh%nnodes))

      ! Delaunay triangulation
      call dtris2 ( green_mesh%nnodes, green_mesh%nodes_coor, &
  &                 green_mesh%nelem, green_mesh%corners, triangle_neighbor )

      ! To know total number of triangular elements
      write(*,*) ' Number of triangular elements is: ', green_mesh%nelem

      ! Write the triangulation to a file=out/elements.out
      call write_triangles(green_mesh%nelem,green_mesh%corners)

      deallocate(triangle_neighbor)
      endsubroutine read_triang_elem



      !==============================================================!
      ! Subroutine read_time(green_mesh)
      ! Reads the Green functions for the problem
      ! according to parameters inside the green_mesh
      ! structure.
      !==============================================================!

      subroutine read_time(green_mesh)

      use triangles
      use interpolate_grid
       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer bit, iunit, i, divide
       integer*8 reclent, col1, ini, fin

       !size of real variable
       bit = 4

       !Divide input data in packages
       divide = 4
       !Read the traction vectors
       iunit=16
       if ( green_mesh%opt_mesh .eq. 1) then
         ! Total elements used for traction matrix (frequency)
         ! For regular cartesian meshes option=1 or 2
         col1=green_mesh%nsstk*green_mesh%nsdip*green_mesh%stcomp*green_mesh%ncomp/divide
         !record length
         reclent=bit*green_mesh%trac_i*col1

        write(*,*) ' Reading Green functions from dat/TRACT_time.bin '
            OPEN(iunit,file='dat/TRACT_time.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
          !read stress-state tensors by columns
          do i=1,divide
            ini = 1+(col1*(i-1))
            fin = col1*i
            write(*,*) ' Read Green s functions: ',i,'/', divide
            read(iunit,rec=i) green_mesh%tractionvec(:,ini:fin)
          enddo
       !to read interpolated stress-state tensors
       elseif ( green_mesh%opt_mesh .eq. 2) then
         ! Total elements used for traction matrix (frequency)
         ! For regular cartesian meshes option=1 or 2
         col1=green_mesh%nsstk*green_mesh%nsdip*green_mesh%stcomp*green_mesh%ncomp/divide
         !record length
         reclent=bit*green_mesh%trac_i*col1

        write(*,*) ' Reading Green s functions from dat/TRACT_time_interp.bin '
            OPEN(iunit,file='dat/TRACT_time_interp.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
          do i=1,divide
            ini = 1+(col1*(i-1))
            fin = col1*i
            write(*,*) ' Read Green s functions:', i, '/', divide
            read(iunit,rec=i) green_mesh%tractionvec(:,ini:fin)
          enddo
       !if iregular mesh is used
       elseif ( green_mesh%opt_mesh .eq. 3) then
          !Read underlying mesh i=1 and interpolate it
          i=1
          call interpolate_green_ref(green_mesh,i)
          !Determine triangulation of mesh
          call read_triang_elem(green_mesh)
          !Determine gauss points for surface integration
          green_mesh%xw = triangle_gauss_points(green_mesh%n_gauss,green_mesh%siz_gauss)
       endif
       close(iunit)


      end subroutine read_time
