  
  implicit none  
  include 'optim_type.h'  
  
  integer :: n                                ! dimension of the problem
  real :: fcost                               ! cost function value
  real,dimension(:),allocatable :: x          ! current point
  real,dimension(:),allocatable :: grad       ! gradient
  real,dimension(:),allocatable :: grad_preco ! preconditioned gradient
  type(optim_type) :: optim                   ! data structure 
  character*4 :: FLAG                         ! communication FLAG
  
