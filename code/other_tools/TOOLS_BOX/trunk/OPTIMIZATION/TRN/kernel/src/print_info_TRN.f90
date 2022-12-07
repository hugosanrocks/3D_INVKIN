!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for printing the     !
! convergence history in the files iterate_TRN.dat    !
! and iterate_TRN_CG.dat                              !
! The file iterate_TRN.dat is similar to the          !
! convergence history files of the other optimization !
! routines of the toolbox (PSTD,PNLCG,LBFGS,PLBFGS)   !
! The file iterate_TRN_CG.dat contains additional     !
! information on the convergence of the inner         !
! conjugate gradient iteration for the computation of !
! the inexact Newton descent direction.               !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine print_info_TRN(n,optim,fcost,FLAG)

  implicit none
  include 'optim_type.h'  

  !IN
  character*4 :: FLAG
  integer :: n
  real :: fcost
  type(optim_type) :: optim !data structure   

  if(optim%print_flag.eq.1) then
     
     if(FLAG.eq.'INIT') then          
        open(10,file='iterate_TRN.dat')     
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A70)') '                                 TRUNCATED NEWTON ALGORITHM                               '
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(10,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(10,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(10,'(A30,ES10.2)') '     Initial norm_grad is   : ',optim%norm_grad
        write(10,'(A30,I7)')     '     Maximum CG iter        : ',optim%niter_max_CG
        write(10,'(A90)') '******************************************************************************************'
        write(10,'(A10,A10,A10,A10,A12,A12,A8,A10,A9,A9)')&
             '   Niter   ' ,&
             '   fk      ' ,&
             '   ||gk||  ' ,&
             '     fk/f0   ' ,&
             '       alpha   ' ,&
             '        nls     ',&
             '  nit_CG',&
             '    eta     ',&
             '  ngrad ',&
             '  nhess '
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8,ES12.2,I8,I8)') &
             optim%cpt_iter,&
             fcost,&
             optim%norm_grad,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%cpt_iter_CG,&
             optim%eta,&
             optim%nfwd_pb,&
             optim%nhess
        open(21,file='iterate_TRN_CG.dat')     
        write(21,'(A90)') '******************************************************************************************'
        write(21,'(A70)') '                                 TRUNCATED NEWTON ALGORITHM                               '
        write(21,'(A70)') '                                      INNER CG HISTORY                                    '
        write(21,'(A90)') '******************************************************************************************'
        write(21,'(A30,ES10.2)') '     Convergence criterion  : ',optim%conv
        write(21,'(A30,I7)')     '     Niter_max              : ',optim%niter_max
        write(21,'(A30,ES10.2)') '     Initial cost is        : ',optim%f0
        write(21,'(A30,ES10.2)') '     Initial norm_grad is   : ',optim%norm_grad
        write(21,'(A30,I7)')     '     Maximum CG iter        : ',optim%niter_max_CG
        write(21,'(A90)') '******************************************************************************************'
     elseif(FLAG.eq.'CONV') then     
        write(10,'(A70)') '**********************************************************************'
        if(optim%cpt_iter.eq.optim%niter_max) then
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
     elseif(optim%comm.eq.'DESC') then
        if(optim%CG_phase.eq.'INIT') then
           write(21,'(A90)') '-------------------------------------------------------------------------------------------------'
           write(21,'(A20,I4,A10,ES12.2)') ' NONLINEAR ITERATION ',optim%cpt_iter, ' ETA IS : ',optim%eta
           write(21,'(A90)') '-------------------------------------------------------------------------------------------------'
           write(21,'(A12,A10,A10,A20)')&
                '  Iter_CG      ',&
                '  qk           ',&
                ' norm_res     ',&
                ' norm_res/||gk||  '
           write(21,'(I8,ES12.2,ES12.2,ES12.2)') &
                optim%cpt_iter_CG,&
                optim%qk_CG,&
                optim%norm_residual,&
                optim%norm_residual/optim%norm_grad       
        else
           write(21,'(I8,ES12.2,ES12.2,ES12.2)') &
                optim%cpt_iter_CG,&
                optim%qk_CG,&
                optim%norm_residual,&
                optim%norm_residual/optim%norm_grad
        endif
     else
        write(10,'(I6,ES12.2,ES12.2,ES12.2,ES12.2,I8,I8,ES12.2,I8,I8)') &
             optim%cpt_iter,&          
             fcost,&
             optim%norm_grad,&
             fcost/optim%f0,&
             optim%alpha,&
             optim%cpt_ls,&
             optim%cpt_iter_CG,&
             optim%eta,&
             optim%nfwd_pb,&
             optim%nhess
     endif
  endif
end subroutine print_info_TRN

