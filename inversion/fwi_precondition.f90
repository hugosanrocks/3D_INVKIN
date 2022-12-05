       subroutine prepare_preconditioner(green_mesh)

       implicit none
!      Define all the variables needed to read models and
!      associated gradient, variables in ../include/green.h
       INCLUDE 'green.h'
!      optim type strcuture
       TYPE (mesh) :: green_mesh

!      Variables needed only here
       INTEGER :: m, n, n2, i, j, k
       INTEGER :: cols, lenrec, ini, fin, iunit1, iunit2
       REAL,DIMENSION(:,:),ALLOCATABLE :: V, SI, VSI, VSIUT, UT, xsol, ATb
       REAL alpha, beta
       INTEGER option

       n = green_mesh%coarse1d           ! coarse grid modelsize
       option = 2                        ! Compute H^-1 and save it in a file
       iunit1 = 25
       iunit2 = 26

       !Estimate preconditioner H^-1 from eigen vector and values
       if ( option .eq. 1)  then

         allocate(SI(n,n),VSI(n,n),VSIUT(n,n),xsol(n,1))
         allocate(V(n,n),UT(n,n),ATb(n,1))
         SI(:,:)=0.
         V(:,:)=0.
         UT(:,:)=0.
         VSI(:,:)=0.
         VSIUT(:,:)=0.

         !Read singular vectors V and UT
         open(iunit1,file='dat/V.matrix',status='old',action='read')
         open(iunit2,file='dat/UT.matrix',status='old',action='read')
         do i=1,n
           do j=1,n
            read(iunit1,*) V(i,j)
            read(iunit2,*) UT(i,j)
           enddo
          enddo
         close(iunit1)
         close(iunit2)

         open(iunit1,file='dat/singular_values.txt',status='old',&
  &           action='read')
         do j=1,7000
          read(iunit1,*) Si(j,j)
          SI(j,j) = 1. / SI(j,j)
         enddo
         close(iunit1)
         write(*,*) ' Singular values and vectors read'

          alpha = 1.
          beta = 0.
          m = n
          k = n
          n = n
          !Forward modeling     Ax = b
          !V*S^-1
          call sgemm('N','N',m,n,k,alpha,V,m,SI,k,beta,VSI,m)
          !V*S^-1*U^T
          print *, 'mul1'
          m = n
          k = n
          n = n
          call sgemm('N','N',m,n,k,alpha,VSI,m,UT,k,beta,VSIUT,m)
          !V*S^-1*U^T*A^T b
          print *, 'mul2'

          green_mesh%precon = VSIUT

          cols = green_mesh%coarse1d/green_mesh%msub
          lenrec = cols*green_mesh%coarse1d
          open(iunit1,file=green_mesh%dat//'preconditioner.bin',status='unknown',&
  &            action='write',access='direct',form='unformatted',recl=lenrec)
          do i=1,green_mesh%msub
            ini = 1 + (cols*(i-1))
            fin = cols*i
            write(iunit1,rec=i) green_mesh%precon(:,ini:fin)
          enddo
          close(iunit1)

          deallocate(SI,VSI,VSIUT)
          deallocate(UT,V)

        !Read pre-computed preconditioner H^-1
        elseif ( option .eq. 2 ) then

          !Read preconditioner H^-1 from  dat/preconditioner.bin
          cols = green_mesh%coarse1d/green_mesh%msub
          lenrec = cols*green_mesh%coarse1d
          open(iunit1,file=green_mesh%dat//'preconditioner.bin',status='old',&
  &            action='read',access='direct',form='unformatted',recl=lenrec)
          do i=1,green_mesh%msub
            ini = 1 + (cols*(i-1))
            fin = cols*i
            read(iunit1,rec=i) green_mesh%precon(:,ini:fin)
          enddo
          close(iunit1)

        else 
          write(*,*) ' Wrong option for preconditioner, check fwi_precondioner.f90'
        endif


        endsubroutine prepare_preconditioner



       subroutine preconditioner(green_mesh,optim)

       implicit none
!      Define all the variables needed to read models and
!      associated gradient, variables in ../include/green.h
       INCLUDE 'green.h'
!      optim type strcuture
       TYPE (mesh) :: green_mesh
       include 'optim_type.h'
       type (optim_type) :: optim

!      Variables needed only here
       INTEGER :: m, n, n2, i, k, iunit1, iunit2
       REAL,DIMENSION(:,:),ALLOCATABLE :: ATb, xsol
       REAL alpha, beta

       n = green_mesh%coarse1d                        ! coarse grid modelsize

       allocate(xsol(n,1))
       allocate(ATb(n,1))

       !Read singular vectors V and UT
       !ATb(:,:) = 0.
       !open(iunit1,file='dat/atb.bin',status='old',&
!  &    action='read',form='unformatted',access='direct',&
!  &    !recl=n)
!       read(iunit1,rec=1) ATb(:,1)
       !close(iunit1)
       !xsol = ATb
       !ATb = xsol * -1.
       xsol(:,1) = 0.
       ATb(:,1) = green_mesh%q_dw(:)

        alpha = 1.
        beta = 0.
        m = n
        k = n
        n2 = 1
        call sgemm('N','N',m,n2,k,alpha,green_mesh%precon,m,ATb,k,beta,xsol,m)
        !do i=1,n
        ! write(63,*) xsol(i,1)
        !enddo

        green_mesh%q_dw(:) = xsol(:,1)

        deallocate(xsol,ATb)
        endsubroutine preconditioner


       subroutine map_interp(green_mesh,optim)

       !This routine sets the relative spacing coordinates
       !along a 2D grid of nodes. This subroutines creates
       !the arrays needed to perfom the bilinear interpolation
       use interplib
       implicit none
       INCLUDE 'green.h'
!      optim type strcuture
       TYPE (mesh) :: green_mesh
       include 'optim_type.h'
       type (optim_type) :: optim

       integer siz, i, j, k, iunit, nsubf, nt
       integer sizint, nsubfint
       integer inx, iny, nxf, nyf, nxint, nyint
       real t, dt
       real, dimension(:,:), allocatable :: coor, coorint, model, modelint
       real, dimension(:,:), allocatable :: snap, snapint
       real, dimension(:), allocatable :: xf, yf, xint, yint
       real xwant, ywant, value
       real dx, dy, x_i, y_i, xint_i, yint_i, dxint, dyint

       iunit=25
       open(iunit,file=green_mesh%dat//'map_interp.info',status='old',&
  &         action='read')

       read(iunit,*) nsubf, nxf, nyf
       read(iunit,*) dx, dy, x_i, y_i
       read(iunit,*) nsubfint, nxint, nyint
       read(iunit,*) dxint, dyint, xint_i, yint_i
       read(iunit,*) nt, dt
       close(iunit)

       !Sizes of fine and coarse model
       siz = green_mesh%modelsize1
       sizint = nsubfint * nt
       allocate(model(siz,1),modelint(sizint,1))
       model(:,1) = optim%q_plb(:)

       allocate(snap(nxf,nyf),snapint(nxint,nyint))
       allocate(xf(nxf),yf(nyf),xint(nxint),yint(nyint))

       !Build X and Y vectors for 2D interpolation
       do i = 1, nxf
         xf(i) = x_i
         x_i = x_i + dx
       enddo
       do i = 1, nyf
         yf(i) = y_i
         y_i = y_i + dy
       enddo
       print *, xf(1), xf(nxf), yf(1), yf(nyf)

       !Build X and Y vectors for 2D interpolation
       do i = 1, nxint
         xint(i) = xint_i
         xint_i = xint_i + dxint
       enddo
       do i = 1, nyint
         yint(i) = yint_i
         yint_i = yint_i + dyint
       enddo
       print *, xint(1), xint(nxint), yint(1), yint(nyint)

       !Interpolation at each snapshot
       do j = 1 , 1                  !each component 1D or 2D
        do i = 1 , nt
        !Arrange slip as snapshots
        k = i
         do iny = 1, nyf
          do inx = 1,nxf               !each subfault
            snap(inx,iny) = model(k,j)
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
           modelint(k,j) = snapint(inx,iny)
           k = k + nt
          enddo
         enddo
!         write(44,*) snapint(28,12)
         if (i .eq. 40) then
         do iny = 1, nyf
          do inx = 1,nxf               !each subfault
            write(55,*) snap(inx,iny)
          enddo
         enddo
         do iny = 1, nyint
          do inx = 1,nxint               !each subfault
            write(54,*) snapint(inx,iny)
          enddo
         enddo
         endif
         !---------------------------------------------------!
        enddo
       enddo



       deallocate(xf,yf,xint,yint)
       deallocate(snap,snapint)
       deallocate(model,modelint)

       endsubroutine map_interp



       subroutine map_interp_time(green_mesh,optim,up_dw)

       !This routine sets the relative spacing coordinates
       !along a 2D grid of nodes. This subroutines creates
       !the arrays needed to perfom the bilinear interpolation
       use interplib
       implicit none
       INCLUDE 'green.h'
!      optim type strcuture
       TYPE (mesh) :: green_mesh
       include 'optim_type.h'
       type (optim_type) :: optim

       integer,intent(inout) :: up_dw
       integer siz, sizint, nsubf, nsubfint, nxf, nyf, nt, ntint
       integer i, j, k, k2, iunit, inx, iny, swichint
       real t, dt, dtint, swich
       real, dimension(:), allocatable :: slip, slipint
       real, dimension(:), allocatable :: snap, snapint
       real, dimension(:), allocatable :: t_orig, t_int
       real twant, value


       iunit=25
       open(iunit,file=green_mesh%dat//'map_interp.info',status='old',&
  &         action='read')

       read(iunit,*) nsubf, nxf, nyf    !This lines are not used
       read(iunit,*) nsubfint    !
       read(iunit,*) nsubfint    !
       read(iunit,*) nsubfint    !
       read(iunit,*) nt, dt, ntint, dtint
       close(iunit)


       if ( up_dw .eq. 1 ) then      !down sampling

         !Do nothing
         !Read the slip solution
         siz = nsubf*nt
         sizint = nsubf*ntint
         print *, up_dw, 'siz sizint', siz, sizint
         allocate(slip(siz),slipint(sizint))
         slip = optim%q_plb

       elseif ( up_dw .eq. 2 ) then  !up sampling

         !Interchange the dt nt between the time series
         swichint = ntint
         ntint = nt
         nt = swichint
         swich = dtint
         dtint = dt
         dt = swich
         !Read the slip solution
         siz = nsubf*nt
         sizint = nsubf*ntint
         print *, up_dw, 'siz sizint', siz, sizint
         allocate(slip(siz),slipint(sizint))
         slip = green_mesh%q_dw

       elseif ( up_dw .eq. 3 ) then      !down sampling
         !Do nothing
         !Read the slip solution
         siz = nsubf*nt
         sizint = nsubf*ntint
         print *, up_dw, 'siz sizint', siz, sizint
         allocate(slip(siz),slipint(sizint))
         slip = green_mesh%grad1

       elseif ( up_dw .eq. 4 ) then  !up sampling

         !Interchange the dt nt between the time series
         swichint = ntint
         ntint = nt
         nt = swichint
         swich = dtint
         dtint = dt
         dt = swich
         !Read the slip solution
         siz = nsubf*nt
         sizint = nsubf*ntint
         print *, up_dw, 'siz sizint', siz, sizint
         allocate(slip(siz),slipint(sizint))
         slip = green_mesh%q_dw

        endif

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
       print *, t_int(ntint), t_orig(nt), 'time'

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
          do i = 1, ntint
            twant = t_int(i)
            call FIND1D(twant,value,t_orig,snap,nt)
            slipint(k2) = value
            k2 = k2 + 1
          enddo
        enddo
       enddo

       if ( up_dw .eq. 1) then
          green_mesh%q_dw = slipint
       elseif ( up_dw .eq. 2) then
          optim%q_plb = slipint
       elseif ( up_dw .eq. 3) then
          green_mesh%q_dw = slipint
       elseif ( up_dw .eq. 4) then
          green_mesh%q_up = slipint
       endif

       iunit = 88
       !NEXT WRITE COMMANDS USED TO CHECK INTERPOLATION
       !Read the slip solution
       open(iunit,file='dat/q_small.dat',status='unknown')
       do i = 1, sizint
       write(iunit,*) t_int(i), slipint(i)
       enddo
       close(iunit)

       iunit= 89
       open(iunit,file='dat/q_big.dat',status='unknown')
       do i = 1, siz
       write(iunit,*) t_orig(i), slip(i)
       enddo
       close(iunit)

       deallocate(t_orig,t_int)
       deallocate(snap,snapint)
       deallocate(slip,slipint)
       endsubroutine map_interp_time


!       subroutine interpolate_grad(green_mesh)

 !      subroutine map_interp(green_mesh)

       !This routine interpolates the input 3D time-space
       !history (gradient or model) along a finer grid of nodes. 
       !This subroutines creates the arrays needed to perfom 
       !the bilinear interpolation in 2D(space)  and linear in 1D (time).

  !     implicit none
  !     INCLUDE 'green.h'
  !     TYPE (mesh) :: green_mesh

  !     integer siz, i, j, k, iunit, nsubf, nt
  !     integer sizint, nsubfint, ntint
  !     integer inx, iny, nxf, nyf, nxint, nyint
  !     real t, dt, dtint, dx, dxint, dy, dyint
  !     real xini, yini, tini
  !     real, dimension(:,:), allocatable :: snap, snapint
   !    real, dimension(:), allocatable :: xf, yf, xint, yint, vector1, vector2
   !    real xwant, ywant, value

!       iunit=25
!       open(iunit,file=green_mesh%dat//'map_interp.info',status='old',&
!  &         action='read')

!       read(iunit,*) nsubf, nxf, nyf
!       read(iunit,*) dx, dy, x_i, y_i
!       read(iunit,*) nsubfint, nxint, nyint
!       read(iunit,*) dxint, dyint, xint_i, yint_i
!       read(iunit,*) nt, dt
!       close(iunit)




!       endsubroutine interpolate_grad















