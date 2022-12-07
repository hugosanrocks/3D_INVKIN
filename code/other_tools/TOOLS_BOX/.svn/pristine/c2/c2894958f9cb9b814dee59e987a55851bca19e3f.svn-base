!*******************************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX                    *!
!*******************************************************!
! This routine is used to generate formatted            !
! ascii output files containing information             !
! on the convergence history.                           !
! The name of these files depend on the                 !
! optimization routine used:                            !
! - steepest descent                : iterate_ST.dat    !
! - nonlinear conjugate gradient    : iterate_CG.dat    !
! - l-BFGS                          : iterate_LB.dat    !
! - preconditioned l-BFGS           : iterate_PLB.dat   !
!                                                       !
! The truncated Newton method and the preconditioned    ! 
! truncated Newont method use their own printing        !
! subroutines to generate their output files since      !
! their format slightly differ.                         !
!                                                       !
! This printing routine is called basically at the      !
! initialization of the optimization and then at each   !
! iteration of the solver                               !
!-------------------------------------------------------!
! INPUT  : character*2 routine_type (solver)            !
!          character*4 FLAG (phase of the solver:       !
!                            initialization,            !
!                            iteration,etc...           !
!          integer n (dimension of the problem)         !
!          optim_type optim (data structure)            !
!-------------------------------------------------------!
subroutine print_info(n,routine_type,optim,fcost,FLAG)

  implicit none
  include 'optim_type.h'  

  !IN
  character*2 :: routine_type
  character*4 :: FLAG
  integer :: n
  real :: fcost
  type(optim_type) :: optim !data structure   
  !Local variables
  real :: ng

  if(optim%print_flag.eq.1) then     
     call normL2(n,optim%grad,ng)
     if(FLAG.eq.'INIT') then
        if(routine_type.eq.'ST') then
           open(10,file='iterate_ST.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '         STEEEPEST DESCENT ALGORITHM         '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'CG') then
           open(10,file='iterate_CG.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '     NONLINEAR CONJUGATE GRADIENT ALGORITHM  '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'LB') then
           open(10,file='iterate_LB.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '             l-BFGS ALGORITHM                '
           write(10,'(A70)') '**********************************************************************'
        elseif(routine_type.eq.'PL') then
           open(10,file='iterate_PLB.dat')     
           write(10,'(A70)') '**********************************************************************'
           write(10,'(A50)') '             PRECONDITIONED l-BFGS ALGORITHM                '
           write(10,'(A70)') '**********************************************************************'
        endif
        write(10,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(10,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(10,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(10,'(A30,ES10.2)') '     Initial norm_grad is   : ',ng
        write(10,'(A70)') '**********************************************************************'
        write(10,'(A10,A10,A10,A10,A12,A12,A12)')&
             '   Niter   ' ,&
             '   fk      ' ,&
             '   ||gk||  ' ,&
             '     fk/f0   ' ,&
             '       alpha   ' ,&
             '       nls     ',&
             '   ngrad    '
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&          
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb
     elseif(FLAG.eq.'CONV') then
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb
        write(10,'(A70)') '**********************************************************************'
        if(optim%cpt_iter.ge.optim%niter_max) then
           write(10,'(A50)') '  STOP: MAXIMUM NUMBER OF ITERATION REACHED    '
        else
           write(10,'(A50)') '  STOP: CONVERGENCE CRITERION SATISFIED        '
        endif
        write(10,'(A70)') '**********************************************************************'
        close(10)
     elseif(FLAG.eq.'FAIL') then
        write(10,'(A70)') '**********************************************************************'
        write(10,'(A50)') '  STOP: LINESEARCH FAILURE    '
        write(10,'(A70)') '**********************************************************************'
        close(10)     
     else
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             ng,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%nfwd_pb    
     endif
  endif
  
end subroutine print_info

