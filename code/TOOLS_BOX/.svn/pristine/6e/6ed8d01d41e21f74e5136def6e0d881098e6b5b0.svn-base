!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! compute the new pairs of models and gradient for    !
! the inverse Hessian approximation                   !
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         integer cpt_lbfgs counts the number of      ! 
!                           already stored pairs      !
!         integer l maximum number of stored pairs    !
!         real,dimension(n) x current point           !
!         real,dimension(n) x current gradient        !
! OUTPUT : sk,yk updated vectors                      !
!-----------------------------------------------------!
subroutine update_LBFGS(n,cpt_lbfgs,l,x,grad,sk,yk)

  implicit none

  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  !IN/OUT
  real,dimension(n,l) :: sk,yk
  !Local variables
  integer :: i
  
  if(cpt_lbfgs.le.l) then
     !---------------------------------------------------!
     ! if the number of stored pairs does not exceed the !
     ! maximum value, then compute a new pair sk yk and  !
     ! update the counter cpt_lbfgs                      !
     !---------------------------------------------------!
     sk(:,cpt_lbfgs)=x(:)-sk(:,cpt_lbfgs)
     yk(:,cpt_lbfgs)=grad(:)-yk(:,cpt_lbfgs)
     cpt_lbfgs=cpt_lbfgs+1
  else
     !---------------------------------------------------!
     ! otherwise, simply update the lth pair             !
     !---------------------------------------------------!
     sk(:,l)=x(:)-sk(:,l)
     yk(:,l)=grad(:)-yk(:,l)
  endif
  
end subroutine update_LBFGS

