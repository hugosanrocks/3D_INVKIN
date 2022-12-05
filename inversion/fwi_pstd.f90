     subroutine fwi_pstd(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n                                       ! dimension of the problem
      real,dimension(:),allocatable :: grad_preco        ! preconditioned gradient
      character*4 :: FLAG                                ! communication FLAG 

!        Save memory for current model, current gradient, preco gradient
         allocate(grad_preco(green_mesh%modelsize2))

!####### Initialize values ################################        
         n=green_mesh%modelsize2 ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=1e-8         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files

         grad_preco(:)=0.        ! Initial pre-conditioned gradient

         !call read_grad(green_mesh)

!#########################################################
 
         do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))

          call PSTD(n,green_mesh%model2,green_mesh%costa,&
  &            green_mesh%grad2,grad_preco,optim,FLAG)

          !Check if computation of new gradient is needed
          if(FLAG.eq.'GRAD') then
             call new_grad(green_mesh)
             grad_preco(:)=green_mesh%grad2(:)   !Pre-condition gradient = current gradient 
          elseif(FLAG.eq.'NSTE') then  !Continue to next step
             green_mesh%iter=green_mesh%iter+1
             write(green_mesh%iter_i,'(I5.5)') green_mesh%iter
             write(44,*) green_mesh%costa
          endif
        
         enddo

      !Write last slip-rate aproximation
      call write_model(green_mesh)
      call write_syn(green_mesh)

      deallocate(grad_preco)
      end subroutine fwi_pstd





     subroutine fwi_pstd1d(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n
      real,dimension(:),allocatable :: grad_preco    ! current gradient
      character*4 :: FLAG                            ! communication FLAG 

         !Save memory for current model, current gradient, preco gradient
         allocate(grad_preco(green_mesh%modelsize1))

         !Initialize values
         n=green_mesh%modelsize1 ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=5e-8         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20
         optim%bound = 1         ! Turn on the boundaries
         print *, optim%bound, 'BOUND'
         allocate(optim%lb(n),optim%ub(n))
         optim%lb(:) = 0.
         optim%ub(:) = 2.5
         optim%threshold = 0.

         grad_preco(:) = 0.

         !----------------------------------------------------!
         ! optimization loop: while convergence not reached or!
         ! linesearch not failed, iterate                     !
         !----------------------------------------------------!

         do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))

          call PSTD(n,green_mesh%model1,green_mesh%costa,&
  &            green_mesh%grad1,grad_preco,optim,FLAG)

          print *, FLAG
          !Check if computation of new gradient is needed
          if(FLAG.eq.'GRAD') then
             call new_grad(green_mesh)
             grad_preco(:)=green_mesh%grad1(:)   !Pre-condition gradient = current gradient 
          elseif(FLAG.eq.'NSTE') then            !Continue to next step
           
          endif

         enddo


      !Write last slip-rate aproximation
      call write_model(green_mesh)
      call write_syn(green_mesh)

      deallocate(grad_preco)
      deallocate(optim%lb,optim%ub)
      end subroutine fwi_pstd1d

