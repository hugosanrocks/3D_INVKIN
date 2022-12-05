     subroutine fwi_lbfgs(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n                                 ! dimension of the problem
      character*4 :: FLAG                          ! communication FLAG 


         !Model size time samples X 3 spatial components X subfaults
         !given by green_mesh%modelsize2 = along strike and dip

!####### Initialize values ################################        
         n=green_mesh%modelsize2 ! dimension
         FLAG='INIT'             ! first flag
         !optim%niter_max=500      ! maximum iteration number 
         optim%conv=1e-8         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20

         green_mesh%iter=0
         write(green_mesh%iter_i,'(I5.5)') green_mesh%iter

         call read_grad(green_mesh)

!#########################################################

  !----------------------------------------------------!
  ! optimization loop: while convergence not reached or!
  ! linesearch not failed, iterate                     !
  !----------------------------------------------------!
   do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     !call model_check(green_mesh) !new step
     call LBFGS(n,green_mesh%model2,green_mesh%costa,green_mesh%grad2,optim,FLAG)
     if(FLAG.eq.'GRAD') then
        call new_grad(green_mesh)
     endif
     write(55,*) green_mesh%costa
  enddo

      !call model_check(green_mesh)
      !Write last slip-rate aproximation
      call write_model(green_mesh)
      call write_syn(green_mesh)


      end subroutine fwi_lbfgs
