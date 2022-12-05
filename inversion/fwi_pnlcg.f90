     subroutine fwi_pnlcg(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim

!     Variables needed only here
      integer :: n, iunit, iunit2, i, j, jj, k, vecl, ii     ! dimension of the problem
      real,dimension(:,:),allocatable :: y, x            ! current point
      real,dimension(:),allocatable :: grad, grad_preco, model, agrad       ! current gradient
      character*4 :: FLAG                                ! communication FLAG 

      integer :: cpt_iter
      real :: mgrad


!        Save memory for current model, current gradient, preco gradient
         allocate(grad_preco(green_mesh%modelsize2))

!####### Initialize values ################################        
         n=green_mesh%modelsize2                  ! dimension
         FLAG='INIT'             ! first flag
         optim%conv=1e-8         ! tolerance for the stopping criterion
         optim%print_flag=1      ! print info in output files 
         optim%debug=.false.     ! level of details for output files
         optim%l=20

         call read_grad(green_mesh)

         grad_preco(:) = green_mesh%grad2(:)

!#########################################################

  !----------------------------------------------------!
  ! optimization loop: while convergence not reached or!
  ! linesearch not failed, iterate                     !
  !----------------------------------------------------!

   do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))

     call PNLCG(n,green_mesh%model2,green_mesh%costa,green_mesh%grad2,grad_preco,optim,FLAG)
     if(FLAG.eq.'GRAD') then
        call new_grad(green_mesh)
        grad_preco(:) = green_mesh%grad2(:)
     endif

   enddo





      !Write last slip-rate aproximation
      call write_model(green_mesh)
      call write_syn(green_mesh)



      deallocate(grad_preco)
      end subroutine fwi_pnlcg
