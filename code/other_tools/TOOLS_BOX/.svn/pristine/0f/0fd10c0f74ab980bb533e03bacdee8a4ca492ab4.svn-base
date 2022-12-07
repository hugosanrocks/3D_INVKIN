!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is the reverse communication   ! 
! mode linesearch shared by all the           ! 
! optimization routines of the toolbox.       !
!                                             !                            
! This linesearch enforces the Wolfe          !
! conditions: sufficient decrease             !
!             sufficient curvature            !
! The Wolfe conditions can be found in        !
! Nocedal, Numerical Optimization, 2nd        !
! edition, p.33                               !
!                                             !
! The linesearch method implemented here is   !
! based first on a bracketing strategy, then  !
! on a dichotomy algorithm. A full description!
! of this strategy can be found in            !
! Numerical Optimizationn Theoretical and     !
! Practical Aspects, J.F.Bonnans, J.C.Gilbert,!
! C. Lemaréchal, C.A. Sagastizábal,           !
! Springer-Verlag, Universitext               !
!                                             !
!---------------------------------------------!
! INPUT  : integer :: n (dimension)           !
!          real fcost (current cost)          !
!          real,dimension(n) grad             !
! OUTPUT : real,dimension(n) x                !
!          optim_typ optim (data structure)   !
!---------------------------------------------!

subroutine std_linesearch(n,x,fcost,grad,optim)
  implicit none
  
  include 'optim_type.h'  

  !IN
  integer :: n
  real :: fcost
  real,dimension(n) :: grad
  !IN/OUT
  real,dimension(n) :: x
  type(optim_type) :: optim
  !Local variables
  real :: q0,new_alpha
    
  if(optim%first_ls) then
     !---------------------------------------!
     ! FIRST LINESEARCH: initialization step !
     !---------------------------------------!
     optim%fk=fcost
     call scalL2(n,grad,optim%descent,q0)
     optim%q0=q0
     !set the search interval bounds to 0
     optim%alpha_L=0.
     optim%alpha_R=0.          
     optim%task='NEW_GRAD'
     optim%first_ls=.false.
     optim%xk(:)=x(:)
     x(:)=optim%xk(:)+optim%alpha*optim%descent(:)
     !IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(optim%bound.eq.1) then
        call project(n,optim,x)
     endif
     optim%cpt_ls=0
  elseif( (optim%cpt_ls.ge.optim%nls_max) .and. (fcost<optim%fk)) then     
     !-----------------------------------------------------------------------!
     ! if the number of linesearch iteration outreaches the maximum allowed  !
     ! but a decrease of the misfit is produced then accept the steplength   !
     !-----------------------------------------------------------------------!
     optim%task='NEW_STEP'        
     optim%first_ls=.true.
     !Compute new x in the descent direction     
     x(:)=optim%xk(:)+optim%alpha*optim%descent(:)     
     !IF BOUNDS ACTIVATED, PROJECT x(:) INTO TO THE FEASIBLE ENSEMBLE
     if(optim%bound.eq.1) then
        call project(n,optim,x)
     endif
  elseif(optim%cpt_ls.ge.optim%nls_max) then     
     !-----------------------------------------------------------------------!
     ! if the number of linesearch iteration outreaches the maximum allowed  !
     ! without decreasing the misfit then the linesearch has failed          !
     !-----------------------------------------------------------------------!
     optim%task='FAILURE!'
  else
     !-----------------------------------------------------------------------!
     ! If not initialization step and number of linesearch iteration ok      !
     ! then perform one linesearch iteration                                 !
     !-----------------------------------------------------------------------!
     call scalL2(n,grad,optim%descent,optim%q)
     if( (fcost.le.(optim%fk+(optim%m1*optim%alpha*optim%q0)))&
          .and.(optim%q.ge.(optim%m2*optim%q0)) ) then
        !--------------------------------------------------------------------!
        ! First test if the Wolfe conditions are satisfied with              !     
        ! current steplength, if this is the case, linesearch                ! 
        ! ends here with success                                             !
        !--------------------------------------------------------------------!
        optim%task='NEW_STEP'
        optim%first_ls=.true.
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%f0 :',optim%f0
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q :', optim%q
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
     elseif (fcost>(optim%fk+(optim%m1*optim%alpha*optim%q0))) then
        !--------------------------------------------------------------------!
        ! If the first condition is not satisfied then shrink the            !
        ! search interval                                                    !
        !--------------------------------------------------------------------!
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'failure 1'
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
        optim%alpha_R=optim%alpha
        new_alpha=(optim%alpha_L+optim%alpha_R)/2.
        optim%alpha=new_alpha        
        optim%task='NEW_GRAD'
        optim%cpt_ls=optim%cpt_ls+1
     elseif( (fcost.le. (optim%fk+(optim%m1*optim%alpha*optim%q0)))&
          .and.(optim%q<(optim%m2*optim%q0) ) ) then
        !--------------------------------------------------------------------!
        ! If the second condition is not satisfied then shrink the           !
        ! search interval unless the right bound of the search interval      !
        ! as not yet been defined                                            !
        !--------------------------------------------------------------------!
        if(optim%debug) then
           if(optim%print_flag.eq.1) then
              write(10,*) 'failure 2'
              write(10,*) 'fcost :',fcost
              write(10,*) 'optim%fk :',optim%fk
              write(10,*) 'optim%alpha :',optim%alpha
              write(10,*) 'optim%q0 :',optim%q0
              write(10,*) 'optim%q :',optim%q
              write(10,*) 'm1 :',optim%m1
              write(10,*) 'm2 :',optim%m2
              write(10,*) 'cpt_ls is : ',optim%cpt_ls
           endif
        endif
        optim%alpha_L=optim%alpha
        if(optim%alpha_R.ne.0.) then
           new_alpha=(optim%alpha_L+optim%alpha_R)/2.
        else
           new_alpha=optim%mult_factor*optim%alpha
        endif
        optim%alpha=new_alpha        
        optim%task='NEW_GRAD'
        optim%cpt_ls=optim%cpt_ls+1                        
     endif
     !Compute new x in the descent direction
     x(:)=optim%xk(:)+optim%alpha*optim%descent(:)     
     !IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(optim%bound.eq.1) then
        call project(n,optim,x)
     endif
  endif
end subroutine std_linesearch



