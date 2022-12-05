!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for the computation  !
! of the Newton descent direction using a matrix-free ! 
! conjugate gradient solver.                          !
! The algorithm is taken from                         !
! Numerical Optimization, Nocedal, 2nd edition, 2006  !
! Algorithm 5.2 p.112                                 !
! The stopping criterion is taken from                !
! Choosing the forcing term in an Inexact Newton      !
! method, Eisenstat and Walker, 1994,                 !
! SIAM Journal on Scientific Computing 17 (1), 16-32  !
! The routine is implemented in a reverse             !
! communication format. The user is requested to give ! 
! Hessian-vector products at each iteration of the    !
! process through the communication flag              !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_TRN(n,grad,optim,FLAG)
  
  implicit none
  include 'optim_type.h'
  
  !IN
  character*4 :: FLAG                        !communication FLAG
  integer :: n                               !dimension of the problem
  real,dimension(n) :: grad                  !gradient at the current point
  !IN/OUT
  type(optim_type) :: optim                  !data structure
  !Local variables
  real :: dHd,norm_residual_prev,alpha,&
       beta,grad_term,descent_scal_Hd
  real,dimension(:),allocatable :: mgrad
    
  if(optim%CG_phase.eq.'INIT') then
     !-----------------------------------------------------!
     ! if optim%CG_phase is INIT, initialize the conjugate !
     ! gradient process                                    !
     !-----------------------------------------------------!
     optim%residual(:)=grad(:)
     optim%d(:)=-1.*optim%residual(:)
     optim%Hd(:)=0.
     optim%descent(:)=0.     
     optim%qk_CG=0.
     optim%hessian_term=0.     
     call normL2(n,optim%residual,optim%norm_residual)
     optim%conv_CG=.false.    
     optim%cpt_iter_CG=0
     call print_info_TRN(n,optim,0e0,FLAG)     
     optim%CG_phase='IRUN'
  else                         
     !-----------------------------------------------------!
     ! else perform one conjugate gradient iteration       !
     !-----------------------------------------------------!     
     call scalL2(n,optim%d,optim%Hd,dHd)        
     if(dHd<0.) then 
        !-----------------------------------------------------!
        ! if dHd < 0, detection of a negative eigenvalue of   !
        ! the Hessian operator, stop the process              !        
        !-----------------------------------------------------!     
        optim%conv_CG=.true.
        write(21,*) 'Negative curvature'
        if(optim%cpt_iter_CG.eq.0) then                     
           !-----------------------------------------------------!
           ! if this is the first iteration, thenreturn the      !
           ! opposite of the gradient as descent direction       !
           ! (steepest descent direction)                        !
           !-----------------------------------------------------!     
           optim%descent(:)=optim%d(:)            
           !-----------------------------------------------------!
           ! if the debug option is activated, compute the       !
           ! quadratic function minimized during the conjugate   !
           ! gradient process (check is this value decresae      !
           ! throughout the CG iterations )                      !
           !-----------------------------------------------------!     
           if(optim%debug) then 
              allocate(mgrad(n))
              mgrad(:)=-1.*grad(:)              
              call normL2(n,optim%residual,optim%norm_residual)
              alpha=(optim%norm_residual**2)/dHd         
              optim%qkm1_CG=optim%qk_CG
              call scalL2(n,optim%descent,mgrad,grad_term)           
              optim%hessian_term=optim%hessian_term+(alpha**2)*dHd
              optim%qk_CG=-grad_term+0.5*optim%hessian_term  
              deallocate(mgrad)
           endif
        endif
     else 
        !-----------------------------------------------------!
        ! if dHd > 0, then perform one conjugate gradient     !
        ! iteration                                           !
        !-----------------------------------------------------!     
        !Update descent direction
        call normL2(n,optim%residual,optim%norm_residual)
        alpha=(optim%norm_residual**2)/dHd
        optim%descent_prev(:)=optim%descent(:)
        optim%descent(:)=optim%descent(:)+alpha*optim%d(:)
        optim%residual(:)=optim%residual(:)+alpha*optim%Hd(:)  
        !Update CG direction
        norm_residual_prev=optim%norm_residual
        call normL2(n,optim%residual,optim%norm_residual)
        beta=(optim%norm_residual**2)/(norm_residual_prev**2)
        optim%d(:)=-1.*optim%residual(:)+beta*optim%d(:)                
        !Update iteration counter 
        optim%cpt_iter_CG=optim%cpt_iter_CG+1
        !-----------------------------------------------------!
        ! if the debug option is activated, compute the       !
        ! quadratic function minimized during the conjugate   !
        ! gradient process (check is this value decresae      !
        ! throughout the CG iterations )                      !
        !-----------------------------------------------------!        
        if(optim%debug) then
           allocate(mgrad(n))
           mgrad(:)=-1.*grad(:)              
           optim%qkm1_CG=optim%qk_CG
           call scalL2(n,optim%descent,mgrad,grad_term)
           call scalL2(n,optim%descent_prev,optim%Hd,descent_scal_Hd)
           optim%hessian_term=optim%hessian_term+(alpha**2)*dHd+&
                2.*alpha*descent_scal_Hd
           optim%qk_CG=-grad_term+0.5*optim%hessian_term
           deallocate(mgrad)
        endif
        !-----------------------------------------------------!
        ! Check if the Eisenstat stopping critertion is       !
        ! satisfied                                           !
        !-----------------------------------------------------!        
        optim%conv_CG=(&
             (optim%norm_residual.le.(optim%eta*optim%norm_grad)).or.&
             (optim%cpt_iter_CG.ge.optim%niter_max_CG))        
        !-----------------------------------------------------!
        ! Print information on the current CG iteration       !
        !-----------------------------------------------------!
        call print_info_TRN(n,optim,0e0,FLAG)     
     endif
  endif
  
end subroutine descent_TRN

