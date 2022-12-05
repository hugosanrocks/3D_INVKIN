       subroutine slip_time_interp(green_mesh)

       !--------------------------------------------!
       !Subroutine devoted to the interpolation
       !of the slip-rate time histories in time.
       !The goal is to interpolate in time the 
       !time histories to have more samples used
       !to compute the forward problem compared to
       !real number of inverted parameters.
       !--------------------------------------------!

       !COMMON VARIABLES
       use interplib
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer siz, i, ii, j, k, k2, iunit, nsubf, nt, ntint
       integer sizint, nsubfint
       integer inx, i_nsub, nxf, nyf, nxint, nyint
       real t, dt, dtint
       real, dimension(:), allocatable :: slip, slipint
       real, dimension(:), allocatable :: modelint
       real, dimension(:), allocatable :: t_orig, t_int
       real twant, value
       real ser(8), vals(8), p

       nsubf = green_mesh%msub
       nxf = green_mesh%nsdip
       nyf = green_mesh%nsstk
       nt = green_mesh%interp_i
       dt = green_mesh%slipdt
       ntint = green_mesh%ntint
       dtint = green_mesh%intdt
       nsubfint = nsubf
       nxint = nxf
       nyint = nyf

       !Size of interpolated model slip-rate
       sizint = nsubfint * ntint * green_mesh%ncomp

       !Allocate memory for each subfault
       allocate(slip(nt),slipint(ntint))
       allocate(t_orig(nt),t_int(ntint))
       allocate(modelint(sizint))

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
       k = nt
       k2 = 1
       do i_nsub = 1, nsubf
        do i = 1, green_mesh%ncomp
         slip(1:nt) = green_mesh%model(j:k)
         !=========================================!
         !1D interpolation                         !
         !=========================================!
         do ii = 1, ntint
          twant = t_int(ii)
          call FIND1D(twant,value,t_orig,slip,nt)
          green_mesh%modelint(k2) = value
          k2 = k2 + 1
         enddo
         j = j + nt
         k = k + nt
        enddo
       enddo

       iunit = 35
       !Print to check time interpolation
       !open(iunit,file='fo.35',status='unknown')
       !do i = 1, green_mesh%modelsize
       !write(iunit,*) green_mesh%model(i)
       !enddo
       !close(iunit)

       iunit = 36
       !Print to check time interpolation
       !open(iunit,file='fo.36',status='unknown')
       !do i = 1, sizint
       !write(iunit,*) modelint(i)
       !enddo
       !close(iunit)


       deallocate(t_orig,t_int)
       deallocate(modelint)
       deallocate(slip,slipint)
       endsubroutine slip_time_interp
