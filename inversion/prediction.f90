     subroutine prediction(green_mesh)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

      !Variables needed only here
      integer :: n, iunit, i, j, k, l, cont                      ! dimension of the problem
      integer :: predsize, subfaults
      real :: timewin, percent
      real,dimension(:),allocatable :: serie, model, grad   ! current gradient
      character*4 :: FLAG                                        ! communication FLAG 


        if (green_mesh%synwin .ge. green_mesh%wintochange) then
         call identify_past(green_mesh)
         !call detect_motion(green_mesh)
        endif 

        allocate(serie(green_mesh%msub))

        iunit=25
        open(iunit,file='dat/times_rup.dat',status='unknown')
        do i=1,green_mesh%msub
         read(iunit,*) serie(i)
        enddo
        close(iunit)

        green_mesh%idsub(:) = 0              !Order of subfaults to predict

        !time at the source location
        timewin = (real(green_mesh%samp4win(green_mesh%synwin))-1.)*green_mesh%slipdt+0.5 !jap


        print *, ' Rupture time for this window:', timewin
        !samples of time history at the source
        !k = int( timewin / green_mesh%slipdt ) + 1 + 2 !2 sample extra
        !jump to the corresponding snapshot of weighting matrix t=t(k)
        !l = green_mesh%msub * ( k - 1 ) + 1
        cont = 0
        k = 1
        do i=1,green_mesh%nsdip
         do j=1,green_mesh%nsstk
          !if ( green_mesh%rup(l) .lt. 4. ) then    !band follow weight
          if ( serie(k) .lt. timewin ) then
            green_mesh%idsub(k) = k
            cont = cont + 1
          else
          endif
          k = k + 1
          !l = l + 1
         enddo
        enddo
        write(*,*) ' Subfaults used:', cont
        print *, green_mesh%dowin, 'win'
        subfaults = cont

        !if (green_mesh%synwin .ge. 5) then
        ! call prior_future(green_mesh,green_mesh%idsub,serie,cont,timewin)
        !endif

        allocate(green_mesh%nodes(cont))
        k = 1
        cont = 0
        do i=1,green_mesh%nsdip
         do j=1,green_mesh%nsstk
          if ( green_mesh%idsub(k) .ne. 0 ) then
           cont = cont + 1
           green_mesh%nodes(cont) = green_mesh%idsub(k)
           !print *, green_mesh%nodes(cont), k, cont
          endif
          k = k+1
         enddo
        enddo


!===================================================!
!     Forward and adjoint using the full model
!===================================================!

      !call rake_uni(green_mesh)

      !Estimate total misfit of the current time window from 0 syn
      green_mesh%syn(:,:) = 0.d0
      call read_obs(green_mesh)
      call residual(green_mesh)
      green_mesh%costdini = green_mesh%costd
      print *, ' misfit inside used window:', green_mesh%costd

      call forward(green_mesh)

      call adjoint(green_mesh)
      print *, ' misfit inside window after forward prediction:', green_mesh%costd

      if (green_mesh%rake_opt .eq. 1) then
       call model_pri1d(green_mesh)
       percent = real(green_mesh%synwin-2)*green_mesh%percent_pri !0.02before change !0.04siv !0.02jap
       if (green_mesh%synwin .eq. 2) then
         percent = 0.4!1
       endif
      else
       call model_pri(green_mesh)
       percent = real(green_mesh%synwin-2)*green_mesh%percent_pri!0.045
       if (green_mesh%synwin .eq. 2) then
         percent = 0.1
       endif
      endif
      print *, 'percent', percent
      print *, green_mesh%costa, green_mesh%costm, percent, 'costd, costm, precent'
      green_mesh%lam4 = percent*green_mesh%costa/green_mesh%costm

      !Read gradient either 1D or 2D case
      call read_grad(green_mesh)

      call penalty_terms(green_mesh,optim)

      !Arrange slip rate model in 1D vector
      !removing subfaults not to update
      if (green_mesh%rake_opt .eq. 1) then
        predsize = cont * (green_mesh%interp_i)
        allocate(model(predsize),grad(predsize))
        !Depth preconditioner before extraction
        if (green_mesh%depth_opt .eq. 1) then
         call depth_preco_go(green_mesh)
        else
        endif
        call extract_model(model,grad,predsize,green_mesh)
        call plbfgs_prediction1d(model,grad,predsize,green_mesh,optim)
      elseif (green_mesh%rake_opt .eq. 2) then
        predsize = cont*2*(green_mesh%interp_i)
        allocate(model(predsize),grad(predsize))
        !Depth preconditioner before extraction
        if (green_mesh%depth_opt .eq. 1) then
         call depth_preco_go(green_mesh)
        else
        endif
        call extract_model(model,grad,predsize,green_mesh)
        call plbfgs_prediction(model,grad,predsize,green_mesh,optim)
       endif

      !Insert the slip-rates removed from the progressive inversion
      call insert_model(model,grad,predsize,green_mesh)

      if (green_mesh%depth_opt .eq. 1) then
       call depth_preco_back(green_mesh)
      else
      endif

      !if (green_mesh%synwin .ge. 5) then
      ! call prior_future(green_mesh,green_mesh%nodes,serie,subfaults,timewin)
      !endif

      !Write last slip-rate aproximation
      if (green_mesh%flag_d .eq. 'PRED') then
      call write_model(green_mesh)
      print *, 'wrote model'
      print *, green_mesh%synwin, 'synwin'
      call write_syn(green_mesh)
      else
      endif
!pause

      !call rake_uni(green_mesh)

      deallocate(serie)
      deallocate(green_mesh%nodes)
      deallocate(model,grad)
      end subroutine prediction


     !===========================================================!
     !===========================================================!

     subroutine plbfgs_prediction(model,grad,predsize,green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n, i, ii, j, k, cont                                        ! dimension of the problem
      integer, intent(inout) :: predsize
      real,dimension(:),allocatable :: grad_preco                  ! current gradient
      character*4 :: FLAG                                          ! communication FLAG 
      real, intent(inout) :: model(predsize), grad(predsize)
      real conv

      !Save memory for current model, current gradient, preco gradient
      allocate(grad_preco(predsize))

      !Initialize values ================================================!
         n=predsize              ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=4.5e-03        ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20
         optim%bound=1           ! turn on bound
         conv=1.0e-02

         cont = n / (green_mesh%interp_i*2)   !number of subfaults used

         !Boundaries to limit rake
         allocate(optim%lb(n),optim%ub(n))
         if (green_mesh%depth_opt .eq. 2) then
          k=1
          do i=1,cont
           do ii=1,2
            do j=1,green_mesh%interp_i
             if ( ii .eq. 2) then
              optim%lb(k) = -0.25
              optim%ub(k) =  0.25
             else
              optim%lb(k) = -1.*green_mesh%ub
              optim%ub(k) =  0.
             endif
             k=k+1
            enddo
           enddo
          enddo
         else
          k=1
          do i=1,cont
           do ii=1,2
            do j=1,green_mesh%interp_i
             if ( ii .eq. 2) then
              optim%lb(k) = -0.25/(green_mesh%fault(1,3)**green_mesh%depth_coef)
              optim%ub(k) =  0.25/(green_mesh%fault(1,3)**green_mesh%depth_coef)
             else
              optim%lb(k) = -1.*green_mesh%ub/(green_mesh%fault(1,3)**green_mesh%depth_coef)
              optim%ub(k) =  0./(green_mesh%fault(1,3)**green_mesh%depth_coef)
             endif
             k=k+1
            enddo
           enddo
          enddo
         endif
         print *, optim%lb(14), optim%ub(15), 'bounds'
         print *, optim%lb(16), optim%ub(17), 'bounds'
         optim%threshold = 0.

         !Iterations
         if (green_mesh%flag_d .eq. 'DETC') then
         optim%niter_max=green_mesh%niter_max !30-(green_mesh%dowin-1)*3
         else
         optim%niter_max=green_mesh%niter_max !30-(green_mesh%dowin-1)*3
         endif

         grad_preco(:) = grad(:)

      !==================================================================!
      !----------------------------------------------------!
      ! optimization loop: while convergence not reached or!
      ! linesearch not failed, iterate                     !
      !----------------------------------------------------
      print *, green_mesh%costa, 'initial cost'
      print *, green_mesh%costdini, 'total window misfit'
      do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
       call PLBFGS(n,model,green_mesh%costa,grad,grad_preco,optim,FLAG)

         if(FLAG.eq.'GRAD') then
          call new_gradpredict(model,grad,predsize,green_mesh,optim)

          green_mesh%conv = green_mesh%costd / green_mesh%costdini
           print *, 'conv data', green_mesh%conv
           if (( green_mesh%conv .lt. conv).and.(optim%cpt_iter.gt.green_mesh%niter_max)) then
             FLAG='CONV'
           endif


         elseif (FLAG .eq. 'PREC') then

         endif

      enddo

      ! if(FLAG.eq.'GRAD') then
      !  call new_gradpredict(model,grad,predsize,green_mesh,optim)
      ! endif
      !enddo

      !call rake_uni(green_mesh)

      deallocate(grad_preco,optim%lb,optim%ub)
      end subroutine plbfgs_prediction

     !===========================================================!
     !===========================================================!

     subroutine plbfgs_prediction1d(model,grad,predsize,green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n, i, j, k                                        ! dimension of the problem
      integer, intent(inout) :: predsize
      real,dimension(:),allocatable :: grad_preco                  ! current gradient
      character*4 :: FLAG                                          ! communication FLAG 
      real, intent(inout) :: model(predsize), grad(predsize)
      real conv, inicost, subfaults

      !Save memory for current model, current gradient, preco gradient
      allocate(grad_preco(predsize))

      !Initialize values ================================================!
         n=predsize              ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=4.5e-03        ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20
         conv=1.0e-2

         optim%bound = 1
         allocate(optim%lb(n),optim%ub(n))
         optim%lb(:) = 0.
         subfaults = predsize / green_mesh%interp_i
         if (green_mesh%depth_opt .eq. 1) then
          k=1
          do i=1,subfaults
           do j=1,green_mesh%interp_i
            optim%ub(k) = green_mesh%ub / &
  &                       ((green_mesh%fault(green_mesh%nodes(i),3)&
  &                       **green_mesh%depth_coef))
            k=k+1
           enddo
          enddo
         else
           optim%ub(:) = green_mesh%ub
         endif
         print *, green_mesh%ub, 'ub'
         !print *, optim%lb(100), optim%ub(101), 'bounds'
         optim%threshold = 0.

         if (green_mesh%flag_d .eq. 'DETC') then
          optim%niter_max=green_mesh%niter_max !-(green_mesh%dowin-1)*3
         else
          optim%niter_max=green_mesh%niter_max !-(green_mesh%dowin-1)*3
         endif

         grad_preco(:) = grad(:)
      !==================================================================!

      !----------------------------------------------------!
      ! optimization loop: while convergence not reached or!
      ! linesearch not failed, iterate                     !
      !----------------------------------------------------
      print *, green_mesh%costa, 'initial cost'
      print *, green_mesh%costdini, 'total window misfit'
      do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
       call PLBFGS(n,model,green_mesh%costa,grad,grad_preco,optim,FLAG)

         if (FLAG.eq.'GRAD') then
          call new_gradpredict(model,grad,predsize,green_mesh,optim)

          green_mesh%conv = green_mesh%costd / green_mesh%costdini
           print *, 'conv data', green_mesh%conv
           if (( green_mesh%conv .lt. conv).and.(optim%cpt_iter.gt.green_mesh%niter_max)) then
             FLAG='CONV'
           endif

         elseif (FLAG .eq. 'PREC') then

         endif

      enddo


       !if(FLAG.eq.'GRAD') then
       ! call new_gradpredict(model,grad,predsize,green_mesh,optim)
       !endif
      !enddo

      deallocate(grad_preco,optim%lb,optim%ub)
      end subroutine plbfgs_prediction1d


     !===========================================================!
     !===========================================================!

      subroutine new_gradpredict(model,grad,predsize,green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim


!     Variables needed only here
      integer :: n, i, j, k, l, cont, m, iunit     ! dimension of the problem
      integer, intent(inout) :: predsize
      real,dimension(:),allocatable :: g
      real lambda
      real, intent(inout) :: model(predsize), grad(predsize)

      call insert_model(model,grad,predsize,green_mesh)

      !call rake_uni(green_mesh)

      if (green_mesh%depth_opt .eq. 1) then
        call depth_preco_back(green_mesh)
      else
      endif

      call model_norm(green_mesh,optim%cpt_iter)

      !call rake_uni(green_mesh)


      !if (green_mesh%depth_opt .eq. 1) then
      !  call depth_preco_back(green_mesh)
      !else
      !endif

      !compose the slip for convolutions (x,y,z)!
      !call model_c(green_mesh%model2,green_mesh%model,green_mesh%interp_i,green_mesh%msub,green_mesh%slipm)
      call build_model(green_mesh)


      call forward(green_mesh)
      !Compute residuals for all stations and components
      call residual(green_mesh)
      iunit=38
      open(iunit,file=green_mesh%out//'data_norm.out',status='old',access='append')
      write(iunit,*) green_mesh%costa, optim%cpt_iter
      close(iunit)

      !TIME DOMAIN GRADIENT COMPUTATION
      call conadjtime(green_mesh)

      !Read gradient either 2D or 1D
      call read_grad(green_mesh)


      call penalty_terms(green_mesh,optim)

      if (green_mesh%depth_opt .eq. 1) then
       call depth_preco_go(green_mesh)
      else
      endif


      call extract_model(model,grad,predsize,green_mesh)


      end subroutine new_gradpredict



     !===========================================================!
     !===========================================================!


      subroutine insert_model(model,grad,predsize,green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: n, i, j, k, l, cont, m 
      integer, intent(inout) :: predsize
      real, intent(inout) :: model(predsize), grad(predsize)


      !Arrange slip rate model in 1D vector
      !Insert back the subfaults that were removed
      if (green_mesh%rake_opt .eq. 1) then
        l = 1
        do i=1,green_mesh%msub
          cont=1+(i-1)*green_mesh%interp_i
          if ( green_mesh%idsub(i) .ne. 0 ) then
           !print *, 'back mod', cont, 'cut mod id', l, 'sub', i
           !Decompose slip vector along stk and along dip = 2 directions
           do k=1,green_mesh%interp_i
             green_mesh%model1(cont) = model(l)
             green_mesh%grad1(cont) = grad(l)
             cont = cont + 1
             l = l + 1
           enddo
          else
           !print *, 'back mod', cont, 'cut mod', l, 'sub', i
           do k=1,green_mesh%interp_i
             green_mesh%model1(cont) = green_mesh%model1(cont)
             green_mesh%grad1(cont) = green_mesh%grad1(cont)
             cont = cont + 1
           enddo
          endif
        enddo
      elseif (green_mesh%rake_opt .eq. 2) then
        l = 1
        do i=1,green_mesh%msub
          cont=1+(i-1)*2*green_mesh%interp_i
          if ( green_mesh%idsub(i) .ne. 0 ) then
           !print *, 'back mod 1',cont, 'cut mod',l, i
           !Decompose slip vector along stk and along dip = 2 directions
           do j=1,2
            do k=1,green_mesh%interp_i
             green_mesh%model2(cont) = model(l)
             green_mesh%grad2(cont) = grad(l)
             cont = cont + 1
             l = l + 1
            enddo
           enddo
          else
           cont = 1+(i-1)*2*green_mesh%interp_i
           !Decompose slip vector along stk and along dip = 2 directions
           do j=1,2
            do k=1,green_mesh%interp_i
             green_mesh%model2(cont) = green_mesh%model2(cont)
             green_mesh%grad2(cont) = green_mesh%grad2(cont)
             cont = cont + 1
            enddo
           enddo
          endif
        enddo
      endif

      end subroutine insert_model



      subroutine extract_model(model,grad,predsize,green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: n, i, j, k, l, cont, m 
      integer, intent(inout) :: predsize
      real, intent(inout) :: model(predsize), grad(predsize)


      !removing subfaults not to update
      if ( green_mesh%rake_opt .eq. 1 ) then
        l = 1
        do i=1,green_mesh%msub
          if ( green_mesh%idsub(i) .ne. 0 ) then
            cont=1+(i-1)*green_mesh%interp_i
            !print *, i, cont, l
            !Remove subfaults not to be updated
               do k=1,green_mesh%interp_i
                 model(l) = green_mesh%model1(cont)
                 grad(l) =  green_mesh%grad1(cont)
                 cont = cont + 1
                 l = l + 1
              enddo
          else
          endif
        enddo
      elseif ( green_mesh%rake_opt .eq. 2 ) then
        l = 1
        do i=1,green_mesh%msub
          if ( green_mesh%idsub(i) .ne. 0 ) then
            cont=1+(i-1)*2*green_mesh%interp_i
            !print *, i, cont, l
            !Remove subfaults not to be updated
            do j=1,2
               do k=1,green_mesh%interp_i
                 model(l) = green_mesh%model2(cont)
                 grad(l) =  green_mesh%grad2(cont)
                 cont = cont + 1
                 l = l + 1
              enddo
            enddo
          else
          endif
        enddo
      endif

      end subroutine extract_model

