  n=2                     ! dimension
  FLAG='INIT'             ! first flag
  optim%niter_max=10000   ! maximum iteration number 
  optim%conv=1e-8         ! tolerance for the stopping criterion
  optim%debug=.false.     ! level of details for output files
  optim%l=20              ! maximum number of stored pairs used for
                          ! the l-BFGS approximation
  optim%pring_flag=1      ! print information
  optim%bound=1           ! activation of bound constraints
  allocate(optim%ub(n))   ! allocate upper bound
  allocate(optim%lb(n))   ! allocate lower bound
  ub(:)=40.               ! set upper bound
  lb(:)=-40.              ! set lower bound
  optim%threshold=1e-2    ! set tolerance for bound constraints
  
  allocate(x(n),grad(n),grad_preco(n)) ! allocation
  x(1)=1.5                             ! set initial point
  x(2)=1.5                             
  
  call Rosenbrock(x,fcost,grad)        ! initial cost and gradient

