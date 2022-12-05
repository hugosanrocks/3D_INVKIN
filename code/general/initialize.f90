!--------------------------------------------------------------------------!
!            Subroutine initialize                                      
!
!    Allocation of arrays which length does not change regarding the 
!    type of inversion used.
!--------------------------------------------------------------------------!
!    Hugo S. Sanchez Reyes 20/05/2018
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, J. Virieux
!---------------------------------------------------------------------------!

      subroutine initialize(green_mesh)

!----------------------Definition of variables--------------------------------------!
      !COMMON VARIABLES
      implicit none
      !Define all the variables needed to read stresses
      ! and calculate the tractions associated "Green's
      ! functions", variables in include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !local variable
      integer :: iunit     
      !Change debug = .false. not to write output files
      green_mesh%debug = .false. 

      !unit number for files
      iunit = 22


      !number of stations times components
      green_mesh%stcomp=green_mesh%nsta*green_mesh%ncomp

      !if iregular mesh is used opt_mesh = 3
      if ( green_mesh%opt_mesh .lt. 3 ) then
         !nothing to be done
      elseif ( green_mesh%opt_mesh .eq. 3 ) then
         !msub (nodes) = number of nodes
         green_mesh%msub = green_mesh%nnodes
      endif


      !1 Fault coordinates (x,y,z)
      allocate(green_mesh%fault(green_mesh%msub,3))
      call monitor_alloc('1 %fault  ',' allo')
      !2 Distances between nodes
      allocate(green_mesh%distcorr(green_mesh%msub,green_mesh%msub))
      call monitor_alloc('2 %distcoo',' allo')
      !3 array for observed seismograms
      allocate(green_mesh%obs(green_mesh%lenobs,green_mesh%stcomp))
      call monitor_alloc('3 %obs    ',' allo')
      !4 matrix to store weights for recordings
      allocate(green_mesh%cd(green_mesh%stcomp,green_mesh%stcomp))
      call monitor_alloc('4 %cd     ',' allo')
      !5 matrix penalizing edges
      allocate(green_mesh%ce(green_mesh%msub,green_mesh%msub))
      call monitor_alloc('5 %ce     ',' allo')
      !6 matrices needed for penalizing rupture time
      allocate(green_mesh%ct(green_mesh%slipsam,green_mesh%slipsam))
      call monitor_alloc('6 %ct     ',' allo')
      !7 vector for expected rupture times
      allocate(green_mesh%rtimes(green_mesh%msub))
      call monitor_alloc('7 %rtimes ',' allo')
      !8 vector for time samples of rupture times
      allocate(green_mesh%rsamp(green_mesh%msub))
      call monitor_alloc('8 %rsamp  ',' allo')

      !9 1D array used for filtering observations if required
      allocate(green_mesh%tseries(green_mesh%lenobs))
      call monitor_alloc('9 %tseries',' allo')

      !10 Array for inverting different time-space windows if required
      allocate(green_mesh%samwin(green_mesh%stcomp,green_mesh%wininv))
      call monitor_alloc('10%samwin ',' allo')
      !11 array for time samples at each station for a given instant
      allocate(green_mesh%synsam(green_mesh%nsta))
      call monitor_alloc('11%synsam ',' allo')
      !12 array for node ID and time windows
      allocate(green_mesh%idsub(green_mesh%msub))
      call monitor_alloc('12%idsub  ',' allo')
      !13 time window at each node
      allocate(green_mesh%win(green_mesh%msub,12))
      call monitor_alloc('13%win    ',' allo')

      !14 matrix used for 2D laplacian filter
      allocate(green_mesh%la(green_mesh%msub,green_mesh%msub)) 
      call monitor_alloc('14%la     ',' allo')

      !15
      !LAPLACIAN FILTER IN SPACE 2D
      !correlation distances lx, lz, lt
      !local angles phi, theta
      allocate(green_mesh%vector_in(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%vector_out(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%lx(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%lz(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%lt(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%phi(green_mesh%nsdip,green_mesh%nsstk))
      allocate(green_mesh%theta(green_mesh%nsdip,green_mesh%nsstk))
      call monitor_alloc('15%laplaci',' allo')
    
      !If triangular mesh is used
      if ( green_mesh%opt_mesh .eq. 3 ) then
        ! Tractionvector from reference grid
        ! 16 From reference grid traction vectors
        allocate(green_mesh%tractionvec_ref(green_mesh%trac_i,&
  &              green_mesh%stcomp*green_mesh%ncomp* &
  &              green_mesh%nsstk*green_mesh%nsdip))
        call monitor_alloc('16%tracref',' allo')
        ! 17 For triangular grid traction vectors
        allocate(green_mesh%tractionvec(green_mesh%trac_i,&
  &              green_mesh%stcomp*green_mesh%ncomp* &
  &              green_mesh%nnodes))
        call monitor_alloc('17%tractio',' allo')
      ! 18 corners of triangles
        allocate(green_mesh%corners(3,3*green_mesh%nnodes))
        call monitor_alloc('18%corners',' allo')
        ! 19 node coordinates
        allocate(green_mesh%nodes_coor(2,green_mesh%nnodes))
        call monitor_alloc('19%nodcoor',' allo')
        ! 20 vector of mesh along strike
        allocate(green_mesh%stk_coor(green_mesh%nnodes))
        call monitor_alloc('20%stkcoor',' allo')
        ! 21 vector of mesh along dip
        allocate(green_mesh%dip_coor(green_mesh%nnodes))
        call monitor_alloc('21%dipcoor',' allo')
        ! 22 vector of mesh along strike reference
        allocate(green_mesh%stk_ref(green_mesh%msub))
        call monitor_alloc('22%stkref ',' allo')
        ! 23 vector of mesh along dip reference
        allocate(green_mesh%dip_ref(green_mesh%msub))
        call monitor_alloc('23%dipref ',' allo')
        ! 24 gauss points
        allocate(green_mesh%xw(green_mesh%siz_gauss,3))
        call monitor_alloc('24%xw     ',' allo')
      else
        ! For cartesian regular grid
        ! 25 traction vectors
        allocate(green_mesh%tractionvec(green_mesh%trac_i,&
  &            green_mesh%stcomp*green_mesh%ncomp*green_mesh%msub))
        call monitor_alloc('25%tractio',' allo')
      endif

      ! 26 Initialize residual arrays
      allocate(green_mesh%res(green_mesh%lenobs,green_mesh%stcomp))
      call monitor_alloc('26%res    ',' allo')
      green_mesh%interpadj_i = green_mesh%interp_i


      end subroutine initialize




      subroutine write_debug(green_mesh)

      !COMMON VARIABLES
      implicit none
      ! Define all the GLOBAL variables
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer iunit

      iunit = 25
      open(iunit,file=green_mesh%out//' debug_data.ascii',status='unknown')

      IF (green_mesh%debug_i .eq. 1) then
      GO TO 111
      ELSEIF (green_mesh%debug_i .eq. 3) then
      GO TO 333
      ENDIF

 111   write(iunit,*) green_mesh%rt, green_mesh%ot
       write(iunit,*) green_mesh%nsta, green_mesh%ncomp, green_mesh%msub
       write(iunit,*) green_mesh%simsam, green_mesh%simt, green_mesh%simdt
       write(iunit,*) green_mesh%interp_i
       write(iunit,*) 'Focal mechanism', green_mesh%stk, green_mesh%dip, green_mesh%rak
       write(iunit,*) 'Fault area', green_mesh%lenstk,'X',green_mesh%lendip,'km2'
       write(iunit,*) 'Subfault area', green_mesh%stk_s,'X',green_mesh%dip_s,'km2'
       write(iunit,*) 'No. of subfault', green_mesh%msub
       write(iunit,*) 'Seismic Moment', green_mesh%moment
       write(iunit,*) 'Slip rate samples and slipdt', &
  &                   green_mesh%slipsam, green_mesh%slipdt
       close(iunit)
      GO TO 999

 333  write(iunit,*) 'Normal vector', green_mesh%vnorm
      write(iunit,*) 'Slip vector', green_mesh%vslip
      write(iunit,*) 'Strike vector', green_mesh%vstk
      close(iunit)
      GO TO 999

 999  return
      end subroutine write_debug


      subroutine destroy_arrays(green_mesh)

      !COMMON VARIABLES
      implicit none
      ! Define all the GLOBAL variables
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      deallocate(green_mesh%fault)       !1
      call monitor_alloc('1 %fault  ',' deal')
      deallocate(green_mesh%distcorr)
      call monitor_alloc('2 %distcoo',' deal')
      deallocate(green_mesh%obs)         !12
      call monitor_alloc('3 %obst   ',' deal')
      deallocate(green_mesh%cd)
      call monitor_alloc('4 %cd     ',' deal')
      deallocate(green_mesh%ce)
      call monitor_alloc('5 %ce     ',' deal')
      deallocate(green_mesh%ct)          !13,14,15,16
      call monitor_alloc('6 %ct     ',' deal')
      deallocate(green_mesh%rtimes)      !17
      call monitor_alloc('7 %rtimes ',' deal')
      deallocate(green_mesh%rsamp)       !18
      call monitor_alloc('8 %rsamp  ',' deal')
      deallocate(green_mesh%tseries)
      call monitor_alloc('9 %tseries',' deal')
      deallocate(green_mesh%samwin)
      call monitor_alloc('10%samwin ',' deal')
      deallocate(green_mesh%synsam)
      call monitor_alloc('11%synsam ',' deal')
      deallocate(green_mesh%idsub)
      call monitor_alloc('12%idsub  ',' deal')

      deallocate(green_mesh%win)
      call monitor_alloc('13%win    ',' deal')

      !laplacian matrix
      deallocate(green_mesh%la)
      call monitor_alloc('14%la     ',' deal')

      !Laplacian filter
      deallocate(green_mesh%lx,green_mesh%lz,green_mesh%lt)
      deallocate(green_mesh%phi,green_mesh%theta)
      deallocate(green_mesh%vector_in,green_mesh%vector_out)
      call monitor_alloc('15%laplaci',' deal')


      if ( green_mesh%opt_mesh .eq. 3 ) then
        deallocate(green_mesh%tractionvec) !8
        call monitor_alloc('16%tacref',' deal')
        deallocate(green_mesh%tractionvec) !8
        call monitor_alloc('17%tactio',' deal')

        deallocate(green_mesh%corners)
        call monitor_alloc('18%corners',' deal')
        deallocate(green_mesh%nodes_coor)   
        call monitor_alloc('19%nodcoor',' deal')
        deallocate(green_mesh%stk_coor)
        call monitor_alloc('20%stkcoor',' deal')
        deallocate(green_mesh%dip_coor)
        call monitor_alloc('21%dipcoor',' deal')
        deallocate(green_mesh%stk_ref)
        call monitor_alloc('22%stkref ',' deal')
        deallocate(green_mesh%dip_ref)
        call monitor_alloc('23%dipref ',' deal')
        deallocate(green_mesh%xw)
        call monitor_alloc('24%xw     ',' deal')
        !deallocate(green_mesh%values)
        !call monitor_alloc('25%values ',' deal')
        !deallocate(green_mesh%nodes)
        !call monitor_alloc('26%nodes  ',' deal')

      else
        deallocate(green_mesh%tractionvec) !8
        call monitor_alloc('25%tractio',' deal')
      endif

      deallocate(green_mesh%res)
      call monitor_alloc('26%res    ',' deal')


      !Variables defined according to source limitsS
      call destroywin(green_mesh)

      !deallocate(green_mesh%precon)
      !deallocate(green_mesh%q_dw)


      !Arrays of variable focal mechanism
      deallocate(green_mesh%rake_i)
      call monitor_alloc('54%rake_i ',' deal')
      deallocate(green_mesh%vnorm)
      call monitor_alloc('55%vnorm  ',' deal')
      deallocate(green_mesh%vstk)
      call monitor_alloc('56%vstk   ',' deal')
      deallocate(green_mesh%vdip)
      call monitor_alloc('57%vdip   ',' deal')
      deallocate(green_mesh%vslip)
      call monitor_alloc('58%vslip  ',' deal')
      deallocate(green_mesh%vslip2)
      call monitor_alloc('59%vslip2 ',' deal')
      deallocate(green_mesh%slipm)
      call monitor_alloc('60%slipm  ',' deal')
      deallocate(green_mesh%vnorm_i)
      call monitor_alloc('61%vnorm_i',' deal')
      deallocate(green_mesh%rake_lim)
      call monitor_alloc('62%rake_li',' deal')
      deallocate(green_mesh%stk_i)
      call monitor_alloc('63%stk_i  ',' deal')
      deallocate(green_mesh%dip_i)
      call monitor_alloc('64%dip_i  ',' deal')
      deallocate(green_mesh%samp4win)
      call monitor_alloc('65%samp4wi',' deal')

      deallocate(green_mesh%twin)
      call monitor_alloc('66%twin   ',' deal')


      end subroutine destroy_arrays
