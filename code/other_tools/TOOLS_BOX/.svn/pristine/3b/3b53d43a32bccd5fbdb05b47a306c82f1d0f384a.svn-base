!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! save the pairs of models and gradient for        !
! the inverse Hessian approximation                   !
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         integer cpt_lbfgs counts the number of      ! 
!                           already stored pairs      !
!         integer l maximum number of stored pairs    !
!         real,dimension(n) x current point           !
!         real,dimension(n) x current gradient        !
! OUTPUT : sk,yk updated vector                       !
!-----------------------------------------------------!
subroutine save_LBFGS(n,cpt_lbfgs,l,x,grad,sk,yk)
  
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
     ! maximum value, then save x and grad               !
     !---------------------------------------------------!
     sk(:,cpt_lbfgs)=x(:)
     yk(:,cpt_lbfgs)=grad(:)   
  else
     !---------------------------------------------------!
     ! otherwise, erase the oldest pair and save the     !
     ! new one (shift)                                   !
     !---------------------------------------------------!
     do i=1,l-1
        sk(:,i)=sk(:,i+1)
        yk(:,i)=yk(:,i+1)                
     enddo
     sk(:,l)=x(:)
     yk(:,l)=grad(:)
  endif
  
end subroutine save_LBFGS
