     subroutine fwi_plbfgs(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n, i, j, k, ii, iunit
      real,dimension(:),allocatable :: grad_preco    ! current gradient
      character*4 :: FLAG                            ! communication FLAG 
      real conv

         !Save memory for current model, current gradient, preco gradient
         allocate(grad_preco(green_mesh%modelsize2))

         !Initialize values
         n=green_mesh%modelsize2 ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=2.0e-4         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20
         optim%bound=1           ! turn on bounds

         !Boundaries to limit rake
         allocate(optim%lb(n),optim%ub(n))
         if (green_mesh%depth_opt .eq. 2) then
          k=1
          do i=1,green_mesh%msub
           do ii=1,2
            do j=1,green_mesh%interp_i
             if ( ii .eq. 2) then
              optim%lb(k) =  -0.25
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
          do i=1,green_mesh%msub
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
         print *, optim%lb(45), optim%ub(45), ' Bounds 2D'
         print *, optim%lb(39), optim%ub(39), ' Bounds 2D'
         optim%threshold = 0.


         !Depth precondition
         if (green_mesh%depth_opt .eq. 1) then
          call depth_preco_go(green_mesh)
         else
         endif


         grad_preco(:) = green_mesh%grad2(:)

         !----------------------------------------------------!
         ! optimization loop: while convergence not reached or!
         ! linesearch not failed, iterate                     !
         !----------------------------------------------------!

         !Save cost value
         iunit=38
         open(iunit,file=green_mesh%out//'data_norm.out',status='old',access='append')
         write(iunit,*) green_mesh%costd, optim%cpt_iter
         close(iunit)


         do while ( ( (FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))  )
           call PLBFGS(n,green_mesh%model2,green_mesh%costa,&
  &                    green_mesh%grad2,grad_preco,optim,FLAG)

           if(FLAG.eq.'GRAD') then
           !Depth precondition
             if (green_mesh%depth_opt .eq. 1) then
              call depth_preco_back(green_mesh)
             else
             endif

             !call rake_uni(green_mesh)

             !Check model norm
             call model_norm(green_mesh,optim%cpt_iter)

             call new_grad(green_mesh,optim)
             green_mesh%conv = green_mesh%costd / green_mesh%costdini
             print *, ' Data misfit ', green_mesh%costd
             print *, ' Data convergence ', green_mesh%conv
             if ((green_mesh%conv .lt. conv).and.(optim%cpt_iter .gt.optim%niter_max)) then
              FLAG='CONV'
             endif

              if (green_mesh%depth_opt .eq. 1) then
               call depth_preco_go(green_mesh)
              else
              endif
           elseif (FLAG .eq. 'PREC') then

           endif
         enddo
        if (green_mesh%depth_opt .eq. 1) then
         call depth_preco_back(green_mesh)
        else
        endif

             call model_norm(green_mesh,optim%cpt_iter)


      !Write last slip-rate aproximation
      call write_model(green_mesh,green_mesh%model2,green_mesh%modelsize2)
      call write_syn(green_mesh)

      deallocate(optim%lb,optim%ub)
      deallocate(grad_preco)
      end subroutine fwi_plbfgs




     subroutine fwi_plbfgs1d(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n, i, j, k, iunit
      integer :: up_dw                               ! Up - Down smapling q_plb
      real,dimension(:),allocatable :: grad_preco    ! current gradient
      character*4 :: FLAG                            ! communication FLAG 
      real conv

         !Save memory for current model, current gradient, preco gradient
         allocate(grad_preco(green_mesh%modelsize1))

         !Initialize values
         n=green_mesh%modelsize1 ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=2.0e-7         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20
         optim%bound = 1         ! Turn on the boundaries
         conv=1.0e-2



         print *, optim%bound, 'Option of bounds 1D'
         allocate(optim%lb(n),optim%ub(n))
         optim%lb(:) = green_mesh%lb
         !Depth precondition
         if (green_mesh%depth_opt .eq. 1) then
           optim%ub(:) = green_mesh%ub/(green_mesh%fault(1,3)**green_mesh%depth_coef)
         else
           optim%ub(:) = green_mesh%ub
         endif
         !print *, optim%lb(:), 'l bounds'
         !print *, optim%ub(:), 'u bounds'
         optim%threshold = 0.

! If precondition and interpolation used
         !Only at first iteration
         ! up_dw = 3
         ! call map_interp_time(green_mesh,optim,up_dw)
         ! call preconditioner(green_mesh,optim)
         ! up_dw = 4
         ! call map_interp_time(green_mesh,optim,up_dw)
         ! grad_preco = green_mesh%q_up

! IF interpolation not needed
!          green_mesh%q_dw = green_mesh%grad1
!          do i=1,n
!           write(55,*) green_mesh%q_dw(i)
!          enddo
!          grad_preco = green_mesh%q_dw

          !Depth precondition
          if (green_mesh%depth_opt .eq. 1) then
           call depth_preco_go(green_mesh)
          else
          endif

          ! IF NO PRECONDITION USED
          grad_preco = green_mesh%grad1

         !----------------------------------------------------!
         ! optimization loop: while convergence not reached or!
         ! linesearch not failed, iterate                     !
         !----------------------------------------------------!

         !Save cost value
         iunit=38
         open(iunit,file=green_mesh%out//'data_norm.out',status='old',access='append')
         write(iunit,*) green_mesh%costd, optim%cpt_iter
         close(iunit)



         do while ( ( (FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))  )

          do i=1,green_mesh%modelsize1
           write(55,*) green_mesh%grad1(i)
           write(56,*) green_mesh%model1(i)
          enddo



         call PLBFGS(n,green_mesh%model1,green_mesh%costa,&
  &                    green_mesh%grad1,grad_preco,optim,FLAG)

           if (FLAG .eq. 'GRAD') then
             !Depth precondition
             if (green_mesh%depth_opt .eq. 1) then
              call depth_preco_back(green_mesh)
             else
             endif
             call model_norm(green_mesh,optim%cpt_iter)
             call new_grad(green_mesh,optim)
             green_mesh%conv = green_mesh%costd / green_mesh%costdini
             print *, ' Data misfit', green_mesh%costd
             print *, ' Data convergence', green_mesh%conv
             if ((green_mesh%conv .lt. conv).and.(optim%cpt_iter .gt.optim%niter_max)) then
              FLAG='CONV'
             endif
             !Filter in space
!             call smooth_preco(green_mesh)
              if (green_mesh%depth_opt .eq. 1) then
               call depth_preco_go(green_mesh)
              else
              endif
           elseif (FLAG .eq. 'PREC') then
!             up_dw = 1
!             call map_interp_time(green_mesh,optim,up_dw)
!             call preconditioner(green_mesh,optim)
!             up_dw = 2
!             call map_interp_time(green_mesh,optim,up_dw)

!             green_mesh%q_dw = optim%q_plb

!             optim%q_plb = green_mesh%q_dw
           endif
           open(iunit,file=green_mesh%out//'data_norm.out',status='old',access='append')
           write(iunit,*) green_mesh%costa, optim%cpt_iter
           close(iunit)
           !call model_norm(green_mesh,optim%cpt_iter)
         enddo

        if (green_mesh%depth_opt .eq. 1) then
         call depth_preco_back(green_mesh)
        else
        endif
        call model_norm(green_mesh,optim%cpt_iter)


      !Write last slip-rate aproximation
      call write_model(green_mesh,green_mesh%model2,green_mesh%modelsize2)
      call write_syn(green_mesh)

      deallocate(grad_preco)
      deallocate(optim%lb,optim%ub)
      end subroutine fwi_plbfgs1d



