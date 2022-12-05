       subroutine grad_time_interp(green_mesh)

       !--------------------------------------------!
       !Subroutine devoted to the interpolation
       !of the gradient time histories.
       !The goal is to interpolate in time the 
       !time histories to have less samples used
       !to compute the gradient to have the same
       !number of samples of the inverted parameters.
       !--------------------------------------------!

       !COMMON VARIABLES
       use interplib
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer i, ii, j, k, iunit, nsubf, nt, ntint
       integer nsubfint
       integer i_nsub, nxf, nyf, nxint, nyint
       real t, dt, dtint
       real, dimension(:), allocatable :: grad
       real, dimension(:), allocatable :: t_orig, t_int
       real twant, value

       nsubf = green_mesh%msub
       nxf = green_mesh%nsdip
       nyf = green_mesh%nsstk
       nt = green_mesh%ntint
       dt = green_mesh%intdt
       ntint = green_mesh%interp_i
       dtint = green_mesh%slipdt
       nsubfint = nsubf
       nxint = nxf
       nyint = nyf

       !Size of interpolated gradient
       !Allocate memory for each subfault
       allocate(grad(nt))
       allocate(t_orig(nt),t_int(ntint))

       !Time coordinates needed for interpolation
       t = 0.
       do i=1,nt
         t_orig(i) = t
         t = t + dt
       enddo
       t = 0.
       do i=1,ntint
         t_int(i) = t
         t = t + dtint
       enddo

       !Interpolate each time history at each node
       j = 1
       k = 1
       do i_nsub = 1, nsubf
        do i = 1, green_mesh%ncomp
         grad(:) = green_mesh%tottrac(:,j)
         !=========================================!
         !1D interpolation                         !
         !=========================================!
         do ii = 1, ntint
          twant = t_int(ii)
          call FIND1D(twant,value,t_orig,grad,nt)
          green_mesh%tracint(ii,j) = value
         enddo
         j = j + 1
        enddo
       enddo

       iunit = 35
       !Print to check interpolation
       !open(iunit,file='fo.35',status='unknown')
       !do i = 1, green_mesh%ntint
       !write(iunit,*) green_mesh%tottrac(i,1)
       !enddo
       !close(iunit)

       iunit = 36
       !Print to check interpolation
       !open(iunit,file='fo.36',status='unknown')
       !do i = 1, green_mesh%interp_i
       !write(iunit,*) green_mesh%tracint(i,1)
       !enddo
       !close(iunit)

       deallocate(t_orig,t_int)
       deallocate(grad)
       endsubroutine grad_time_interp
