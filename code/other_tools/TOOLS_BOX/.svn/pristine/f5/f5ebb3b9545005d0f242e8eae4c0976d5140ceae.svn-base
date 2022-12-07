!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! Computation of the descent direction given by the 
! preconditioned nonlinear conjugate gradient algorithm 
! of Dai and Yuan 
! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
! method with a strong global convergence property,   !
! SIAM Journal on Optimization, 10 (1999), pp. 177–182!
!                                                     !
! See also Nocedal, Numerical optimization,           !
! 2nd edition p.132                                   !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real,dimension(n) :: grad,grad_preco       !
! INPUT/OUTPUT : optim_typ optim (data structure)     !
!-----------------------------------------------------!
subroutine descent_PNLCG(n,grad,grad_preco,optim)
  
  implicit none
  include 'optim_type.h'  
  !IN
  integer :: n
  real,dimension(n) :: grad,grad_preco
  !IN/OUT
  type(optim_type) :: optim !data structure   
  !Local variables
  real :: gkpgk,skpk,beta
  real,dimension(:),allocatable :: sk
  
  !------------------------------------------------------------!
  ! Storing old descent direction                              !
  !------------------------------------------------------------!
  optim%descent_prev(:)=optim%descent(:)
  
  !------------------------------------------------------------!
  ! Computation of beta                                        !
  !------------------------------------------------------------!  ! 
  call scalL2(n,grad,grad_preco,gkpgk)  
  allocate(sk(n))
  sk(:)=grad(:)-optim%grad_prev(:)
  call scalL2(n,sk,optim%descent_prev,skpk)
  beta=gkpgk/skpk
  
  !------------------------------------------------------------!
  ! Safeguard (may be useful in some cases)                    !
  !------------------------------------------------------------! 
  if((beta.ge.1e5).or.(beta.le.-1e5)) then     
     beta=0.
  endif
  
  !------------------------------------------------------------!
  ! Computation of the descent direction                       !
  !------------------------------------------------------------! 
  optim%descent(:)=-1.*grad_preco(:)+beta*optim%descent_prev(:)
  
  !------------------------------------------------------------!
  ! Deallocation                                               ! 
  !------------------------------------------------------------!
  deallocate(sk)

end subroutine descent_PNLCG
