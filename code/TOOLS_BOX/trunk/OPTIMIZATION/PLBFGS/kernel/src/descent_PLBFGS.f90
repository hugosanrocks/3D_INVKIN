!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the two parts of the computation !
! of the l-BFGS descent direction following the two   !
! recursion loop algorithm                            ! 
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
! Basically, the two loops have been split, compared  !
! to the standard l-BFGS algorithm, to allow the user !
! to use its preconditioner as an initial estimation  !
! of the inverse Hessian operator                     ! 
!-----------------------------------------------------!

!-----------------------------------------------------!
!FIRST LOOP: descent1_PLBFGS                          !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n,l) sk yk pairs of vectors  !
!                             for l-BFGS approximation!
!         integer cpt_lbfgs current number of stored  !
!                           pairs                     !
!         integer l maximum number of stored pairs    !
! OUTPUT : real,dimension(n) descent                  !
!-----------------------------------------------------!
subroutine descent1_PLBFGS(n,grad,sk,yk,cpt_lbfgs,l,optim)
  
  implicit none
  include 'optim_type.h'

  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  real,dimension(n,l) :: sk,yk
  type(optim_type) :: optim
  !IN/OUT
  real,dimension(n) :: descent
  !Local variables
  real :: beta,gamma,gamma_num,gamma_den
  real :: norml2_yk,norml2_sk
  integer :: i,borne_i
  
  borne_i=cpt_lbfgs-1
  allocate(optim%alpha_plb(cpt_lbfgs))
  allocate(optim%rho_plb(cpt_lbfgs))  
  optim%q_plb(:)=grad(:)
  do i=1,borne_i
     call scalL2(n,yk(:,borne_i-i+1),sk(:,borne_i-i+1),optim%rho_plb(borne_i-i+1))
     optim%rho_plb(borne_i-i+1)=1./optim%rho_plb(borne_i-i+1)
     call scalL2(n,sk(:,borne_i-i+1),optim%q_plb(:),optim%alpha_plb(borne_i-i+1))
     optim%alpha_plb(borne_i-i+1)=optim%rho_plb(borne_i-i+1)*optim%alpha_plb(borne_i-i+1)
     optim%q_plb(:)=optim%q_plb(:)-optim%alpha_plb(borne_i-i+1)*yk(:,borne_i-i+1)
  enddo
  
end subroutine descent1_PLBFGS

!-----------------------------------------------------!
!SECOND LOOP: descent2_PLBFGS                         !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n,l) sk yk pairs of vectors  !
!                             for l-BFGS approximation!
!         integer cpt_lbfgs current number of stored  !
!                           pairs                     !
!         integer l maximum number of stored pairs    !
! OUTPUT : real,dimension(n) descent                  !
!-----------------------------------------------------!
subroutine descent2_PLBFGS(n,grad,sk,yk,cpt_lbfgs,l,optim,descent)
  
  implicit none
  include 'optim_type.h'
  
  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  real,dimension(n,l) :: sk,yk
  type(optim_type) :: optim
  !IN/OUT
  real,dimension(n) :: descent
  !Local variables
  real :: beta,gamma,gamma_num,gamma_den
  real :: norml2_yk,norml2_sk
  integer :: i,borne_i
  
  borne_i=cpt_lbfgs-1  
  call scalL2(n,sk(:,borne_i),yk(:,borne_i),gamma_num)
  call normL2(n,yk(:,borne_i),gamma_den)
  gamma=gamma_num/(gamma_den**2)     
  descent(:)=gamma*optim%q_plb(:)
  do i=1,borne_i
     call scalL2(n,yk(:,i),descent(:),beta)
     beta=optim%rho_plb(i)*beta
     descent(:)=descent(:)+(optim%alpha_plb(i)-beta)*sk(:,i)
  enddo
  descent(:)=-1.*descent(:)
  deallocate(optim%alpha_plb,optim%rho_plb)
  
end subroutine descent2_PLBFGS
