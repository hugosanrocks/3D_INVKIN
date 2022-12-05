       program slip_interp_time


       implicit none
       integer siz, i, j, k, k2, iunit, nsubf, nt, ntint
       integer sizint, nsubfint
       integer inx, iny, nxf, nyf, nxint, nyint
       real t, dt, dtint
       real, dimension(:), allocatable :: slip, slipint
       real, dimension(:), allocatable :: snap, snapint
       real, dimension(:), allocatable :: t_orig, t_int
       real twant, value
       real ser(8), vals(8), p

       !ser = [0., 0.16, 0.32, 0.48, 0.64, 0.8, 0.96, 1.12]
       !vals = ser
       !p = 1.16
       !call FIND1D(p,value,ser,vals,8)
       !print *, p, value

       nsubf = 288
       nxf = 24
       nyf = 12
       nt = 57
       dt = 0.16
       ntint = 225
       dtint = 0.04
       nsubfint = nsubf
       nxint = nxf
       nyint = nyf

       iunit = 26
       !Read the slip solution
       open(iunit,file='dat/xsol.dat',status='old',action='read')
       siz = nsubf * nt
       sizint = nsubfint * ntint

       allocate(slip(siz),slipint(sizint))
       do i=1,siz
        read(iunit,*) slip(i)
       enddo
       close(iunit)

       allocate(snap(nt),snapint(ntint))
       allocate(t_orig(nt),t_int(ntint))

       !Build t_orig and t_int vectors for 1D interpolation
       t = 0.
       do i = 1, nt
         t_orig(i) = t
         t = t + dt
       enddo
       t = 0.
       do i = 1, ntint
         t_int(i) = t
         t = t + dtint
       enddo

       !Interpolate each time history at each node
       k = 1
       k2= 1
       do iny = 1, nyf
        do inx = 1, nxf
          do i = 1 , nt
            snap(i) = slip(k)
            k = k + 1
          enddo
        !===================================================!
        !Interpolate, arrange as snapshot and save as vector
        !if ( iny .eq. 9) then
        ! if (inx .eq. 18) then
          do i = 1, ntint
            twant = t_int(i)
            call FIND1D(twant,value,t_orig,snap,nt)
           slipint(k2) = value
         !  print *, twant, value, i
           k2 = k2 + 1
          enddo
        !endif
       !endif
        enddo
       enddo

       !Read the slip solution
       open(iunit,file='dat/model1dinterp.dat',status='unknown')
       do i = 1, sizint
       write(iunit,*) slipint(i)
       enddo
       close(iunit)



       deallocate(t_orig,t_int)
       deallocate(snap,snapint)
       deallocate(slip,slipint)
       endprogram slip_interp_time
