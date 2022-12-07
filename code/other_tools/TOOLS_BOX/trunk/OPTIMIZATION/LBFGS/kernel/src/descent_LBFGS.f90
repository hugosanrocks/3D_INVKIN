!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routines is used by the l-BFGS routine to      !
! compute l-BFGS descent direction following the two  !
! recursion loop algorithm                            ! 
! Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     !
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n,l) sk yk pairs of vectors  !
!                             for l-BFGS approximation!
!         integer cpt_lbfgs current number of stored  !
!                           pairs                     !
!         integer l maximum number of stored pairs    !
! OUTPUT : real,dimension(n) descent                  !
!-----------------------------------------------------!
subroutine descent_LBFGS(n,grad,sk,yk,cpt_lbfgs,l,descent)
  
  implicit none
  
  !IN
  integer :: n,cpt_lbfgs,l
  real,dimension(n) :: x,grad
  real,dimension(n,l) :: sk,yk
  !IN/OUT
  real,dimension(n) :: descent
  !Local variables
  real,dimension(:),allocatable :: q
  real,dimension(:),allocatable :: alpha,rho
  real :: beta,gamma,gamma_num,gamma_den
  real :: norml2_yk,norml2_sk
  integer :: i,borne_i
  
  borne_i=cpt_lbfgs-1
  
  !------------------------------------------!
  ! SAFEGUARD                                !
  !------------------------------------------!
  call normL2(n,sk(:,borne_i),norml2_sk)
  call normL2(n,yk(:,borne_i),norml2_yk)
  if( (norml2_sk==0.).or.(norml2_yk==0.)) then
     descent(:)=-1.*grad(:)     
  else
     !------------------------------------------!
     ! First phase of the recursion loop        !
     !------------------------------------------!
     allocate(alpha(cpt_lbfgs))
     allocate(rho(cpt_lbfgs))
     allocate(q(n))
     q(:)=grad(:)
     do i=1,borne_i
        call scalL2(n,yk(:,borne_i-i+1),sk(:,borne_i-i+1),rho(borne_i-i+1))
        rho(borne_i-i+1)=1./rho(borne_i-i+1)
        call scalL2(n,sk(:,borne_i-i+1),q(:),alpha(borne_i-i+1))
        alpha(borne_i-i+1)=rho(borne_i-i+1)*alpha(borne_i-i+1)
        q(:)=q(:)-alpha(borne_i-i+1)*yk(:,borne_i-i+1)
     enddo
     call scalL2(n,sk(:,borne_i),yk(:,borne_i),gamma_num)
     call normL2(n,yk(:,borne_i),gamma_den)
     !------------------------------------------!
     ! Scaling by gamma                         !
     !------------------------------------------!
     gamma=gamma_num/(gamma_den**2)     
     descent(:)=gamma*q(:)
     !------------------------------------------!
     ! Second phase of the recursion loop       !
     !------------------------------------------!
     do i=1,borne_i
        call scalL2(n,yk(:,i),descent(:),beta)
        beta=rho(i)*beta
        descent(:)=descent(:)+(alpha(i)-beta)*sk(:,i)
     enddo
     descent(:)=-1.*descent(:)
     deallocate(q,alpha,rho)
  endif
  
end subroutine descent_LBFGS
