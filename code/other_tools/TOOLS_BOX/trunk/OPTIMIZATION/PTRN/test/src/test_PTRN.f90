!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This program provides an example for using the      !
! preconditioned truncated Newton algorithm within the!
! SEISCOPE OPTIMIZATION TOOLBOX.                      !
!                                                     !
! The algorithm implemented is described in           !
! L.Metivier, R. Brossier, J. Virieux, S.Operto,      !
! Truncated Newton and full waveform inversion, 2013  !
! SIAM Journal on Scientific Computing, Vol. 35,      !
! No. 2, pp. B401–B437,                               !         
!                                                     !
! The 2D Rosenbrock function is minimized             !
! f(x,y)  = (1-x)**2+100.*(y-x**2)**2                 !
!                                                     !
! The computation of this function and its gradient   !
! is implemented in                                   !
! ../../../COMMON/test/rosenbrock.f90                 !
! routine rosenbrock                                  !
!                                                     !
! The preconditioned truncated Newton method also     !
! requires to compute Hessian-vector products. This   !
! is implemented in                                   ! 
! ../../../COMMON/test/rosenbrock.f90                 !
! routine rosenbrok_hess                              !
!                                                     !
! The communication with the PTRN optimizer is        ! 
! achieved through the variable FLAG. The only        !
! difference with the TRN optimizer is that a         !
! preconditioner can be used to improve the           !
! convergence speed of the inner conjugate gradient   !
! iterations. This preconditioner is applied to the   !
! gradient, and to the residual of the inner linear   !
! system.  
!                                                     !
! Here are the actions requested for the user         !
! 1. Before first call, the user must set the         !
!    optimization parameters, namely                  !
!    - maximum number of iterations optim%niter_max   !
!    - tolerance for the stopping criterion           !
!      optim%conv                                     !
!      (this stopping criterion is the default one    !
!        in the TOOLBOX and is based on the decrease  !
!        of the misfit function scaled to one at the  !
!        first iteration)                             !
!     - level of details in the output                !
!       * if optim%debug=.true. then the information  !
!         on the linesearch process is provided       !
!         additional information related to the       !
!         quadratic function minimized throughout the !
!         inner conjugate gradient iterations is also !
!         provided as a quality control               !
!       * if optim%debug=.false. a more compact output!
!         file is provided                            !
!   - initialize x with an arbitrary initial guess    !
!   - initialize fcost with the cost function         !
!     associated with the initial guess               !
!   - initialize grad with the gradient associated    !
!     with the initial guess                          !
!   - initialize grad_preco with the gradient         !
!     associated with the initial guess multiplied by !
!     the preconditioner                              !
! 2. On first call, FLAG must be set to INIT          !
! 3. When the FLAG returned by PTRN is 'GRAD', then   !
!    the user must compute the cost, the gradient     !
!    and the preconditioned gradient at current       !
!    iterate x in fcost, grad, grad_preco             !
! 4. When the FLAG returned by PTRN is 'PREC', then   !
!    multiply the vector optim%residual by the        !
!    preconditioner and store the result in           !
!    optim%residual_preco                             !
! 5. When the FLAG returned by PTRN is 'HESS', then   !
!    multiply the vector optim%d by the Hessian       !
!    and store the result in optim%Hd                 ! 
!-----------------------------------------------------!
program test_PTRN
  
  implicit none  
  include 'optim_type.h'  
  
  integer :: n                                ! dimension of the problem
  real :: fcost                               ! cost function value
  real,dimension(:),allocatable :: x          ! current point
  real,dimension(:),allocatable :: x_exact    ! exact solution
  real,dimension(:),allocatable :: grad       ! current gradient
  real,dimension(:),allocatable :: grad_preco ! preconditioned current gradient
  type(optim_type) :: optim                   ! data structure for the optimizer
  character*4 :: FLAG                         ! communication FLAG 

  !----------------------------------------------------------------!
  !variables and parameter for L2-norm
  !----------------------------------------------------------------!
  real :: l2err                  ! L2-norm between x_exact(:)=1e0 and x
  real, parameter :: l2tol=1e-5  ! Tolerance criteria for L2-norm
  

  !----------------------------------------------------!
  ! parameter initialization                           !
  !----------------------------------------------------!
  n=2                     ! dimension
  FLAG='INIT'             ! first flag
  optim%niter_max=100     ! maximum iteration number
  optim%conv=1e-8         ! tolerance for the stopping criterion
  optim%print_flag=1      ! print info in output files 
  optim%debug=.false.     ! level of details for output files
  optim%niter_max_CG=5    ! maximum number of inner conjugate gradient 
                          ! iterations
  
  !----------------------------------------------------!
  ! intial guess                                       !
  !----------------------------------------------------!
  allocate(x(n),grad(n),grad_preco(n))
  x(1)=1.5
  x(2)=1.5
  
  !----------------------------------------------------!
  ! computation of the cost and gradient associated    !
  ! with the initial guess                             !
  !----------------------------------------------------!
  call Rosenbrock(x,fcost,grad)

  !----------------------------------------------------!
  ! multiplication of the gradient by the              !
  ! preconditioner (here preconditioner is identity:   !
  ! no preconditioning)                                !
  !----------------------------------------------------!
  grad_preco(:)=grad(:)  
  
  !----------------------------------------------------!
  ! optimization loop: while convergence not reached or!
  ! linesearch not failed, iterate                     !
  !----------------------------------------------------!
  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     call PTRN(n,x,fcost,grad,grad_preco,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        !----------------------------------------------------!
        ! if FLAG is GRAD, then compute cost, gradient and   !
        ! preconditioned gradient in fcost, grad, grad_preco !
        !----------------------------------------------------!
        call Rosenbrock(x,fcost,grad)        
        grad_preco(:)=grad(:)
     elseif(FLAG.eq.'HESS') then
        !----------------------------------------------------!
        ! if FLAG is HESS, then multiply optim%d by the      !
        ! Hessian operator and store the result in optim%Hd  !
        !----------------------------------------------------!
        call Rosenbrock_Hess(x,optim%d,optim%Hd)
     elseif(FLAG.eq.'PREC') then
        !----------------------------------------------------!
        ! if FLAG is PREC, then multiply optim%residual by   !
        ! preconditioner and store the result in             !
        ! optim%residual_preco                               !
        !----------------------------------------------------!
        optim%residual_preco(:)=optim%residual(:)
     endif
  enddo  

  !L2-Norm
  allocate(x_exact(n))
  x_exact(1)=1e0
  x_exact(2)=1e0

  l2err = 0.
  l2err=sum(x(:)-x_exact(:))**2

  write(*,*) '---------------------------------------------'
  !Helpful console writings
  write(*,*) 'FINAL iterate is : ', x(:)
  write(*,*) 'See the convergence history in :'
  write(*,*) '  iterate_PTRN.dat and iterate_PTRN_CG.dat'
  if (l2err .gt. l2tol) then
    write(*,*) 'CAUTION: Tolerance =   1e-5'
    write(*,*) '         ERR       = ', l2err
    write(*,*) '--- OPTIMIZATION PTRN .............*** Failed'
  else
    write(*,*) '--- OPTIMIZATION PTRN .............*** Passed'
  endif
  write(*,*) '---------------------------------------------'

  deallocate(x_exact)

end program test_PTRN

