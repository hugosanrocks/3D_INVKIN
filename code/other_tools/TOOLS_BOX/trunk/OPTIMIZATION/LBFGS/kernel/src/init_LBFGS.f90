!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is used for the initialization of the  !
! reverse communication mode l-BFGS algorithm from    !
! the TOOLBOX.                                        !
!                                                     !
! The parameters for the linesearch are set here, as  !
! well as the data structure allocations required for !
! the use of the algorithm                            !
!-----------------------------------------------------!
! INPUT : integer n, dimension of the problem         !
!       : real fcost                                  !
!       : real,dimension(n) x,grad first iterate and  !
!                           gradient                  !
! INPUT/OUTPUT : optim_type optim data structure      ! 
!-----------------------------------------------------!
subroutine init_LBFGS(n,l,x,fcost,grad,optim)
  
  implicit none
  include 'optim_type.h'  
  !IN
  integer :: n,l 
  real :: fcost
  real,dimension(n) :: x,grad
  !IN/OUT
  type(optim_type) :: optim !data structure   
  !Local variable
  real,dimension(:),allocatable :: mgrad
  
  !---------------------------------------!
  ! set counters                          !
  !---------------------------------------!
  optim%cpt_iter=0
  optim%f0=fcost
  optim%nfwd_pb=0  
  allocate(optim%sk(n,l),optim%yk(n,l))
  optim%sk(:,:)=0.
  optim%yk(:,:)=0.
  optim%cpt_lbfgs=1
  
  !---------------------------------------!
  ! initialize linesearch parameters      !
  ! by default, the max number of         !
  ! linesearch iteration is set to 20     !
  ! and the initial steplength is set to 1!
  !---------------------------------------! 
  optim%m1=1e-4 ! Wolfe conditions parameter 1 (Nocedal value)
  optim%m2=0.9  ! Wolfe conditions parameter 2 (Nocedal value)
  optim%mult_factor=10 ! Bracketting parameter (Gilbert value)
  optim%fk=fcost
  optim%nls_max=20 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1. ! first value for the linesearch steplength 
  
  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(optim%xk(n))
  optim%xk(:)=x(:)   
  allocate(optim%grad(n))
  optim%grad(:)=grad(:)
  allocate(optim%descent(n))

  !---------------------------------------!
  ! first descent direction               !
  !---------------------------------------!
  optim%descent(:)=-1.*grad(:)  
  call save_LBFGS(n,optim%cpt_lbfgs,l,x,grad,optim%sk,optim%yk)
  
end subroutine init_LBFGS
