!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- OPTIMIZATION OPTION SUBROUTINE -----------
!
!     This program:
!
!    4) Optimization:
!
!     4.1) Searchs for a slip-rate model that minimizes the 
!          cost function (L2 data). Four strategies can be
!          selected: Precondition Steepest descent (PSTD),
!          Limited Memory BFGS (LBFGS), Precondition LBFGS (PLBFGS)
!          and the Precondition Nonlinear Conjugate Gradient (PNLCG).
!
!          dat/fwioption.info
!          1 PSTD
!          2 LBFGS
!          3 PLBFGS
!          4 PNLCG
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------

       subroutine fwi_option(green_mesh,optim)

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      !Optimization TOOL_BOX
      include 'optim_type.h'
      type (optim_type) :: optim

      !Local variables
      integer iunit, fwiopt

      !Use any mehtod to follow the -gradient direction to find
      !a model (solution) that minimizes the cost function

      !Read the option selected
      iunit=12
      open(iunit,file=green_mesh%dat//'fwioption.info',status='unknown')
      read(iunit,*) green_mesh%fwiopt
      read(iunit,*) optim%niter_max
      read(iunit,*) green_mesh%weig
      read(iunit,*) green_mesh%lam1, green_mesh%lam2, green_mesh%lam3
      close(iunit)

      if (green_mesh%fwiopt .eq. 1) then
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PSTD INITIAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
        if (green_mesh%rake_opt .eq. 1) then
          write(*,*) ' Call to pstd1d'
          call fwi_pstd1d(green_mesh,optim)
        elseif (green_mesh%rake_opt .eq. 2) then
          write(*,*) ' Call to pstd'
          call fwi_pstd(green_mesh,optim)
        endif
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PSTD FINAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
      else if(green_mesh%fwiopt .eq. 2) then
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: LBFGS INITIAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
        call fwi_lbfgs(green_mesh,optim)
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: LBFGS FINAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
      else if(green_mesh%fwiopt .eq. 3) then
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PLBFGS INITIAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
        if (green_mesh%rake_opt .eq. 1) then
          write(*,*) ' Call to plbfgs1d'
          call fwi_plbfgs1d(green_mesh,optim)
        elseif (green_mesh%rake_opt .eq. 2) then
          write(*,*) ' Call to plbfgs'
          call fwi_plbfgs(green_mesh,optim)
        endif
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PLBFGS FINAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
      else if(green_mesh%fwiopt .eq. 4) then
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PNLCG INITIAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
        call fwi_pnlcg(green_mesh,optim)
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIMIZTION OPTION: PNLCG FINAL COST ',green_mesh%costa
        WRITE(6, *) '================================================='
      else
        WRITE(6 ,*) ''
        WRITE(6, *) '================================================='
        WRITE(6, *) ' OPTIM OPTION NOT VALID, CHECK dat/fwioption.info '
        WRITE(6, *) '================================================='
      endif

      end subroutine fwi_option
