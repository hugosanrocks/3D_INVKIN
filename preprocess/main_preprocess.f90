!----------------------------------------------------------------------
!
!     Preprocess of stress tensors to compute
!     unitary traction vectors along strike and dip
!     directions of a fault geometry.
!
!     ------------- Procedure info ---------------------------
!
!     This program:
!
!     1) Preprocess:
!
!      1.1) Reads the: 
!              - focal mechanism information.
!              - simulation (Green functions) information.
!              - slip rate model information to reconstruct.
!      1.2) Using this information, the program 
!           performs a coordinate transformation
!           slip rate (x,y,z).
!      1.3) Prepares the stress tensor by:
!             - Reading files P_CX TAU_*_CX SIGMA_**_CX.
!             - Rearranging this files into a unique file
!               for each station and component.
!             - Interpolating at the stress tensor at 
!               sliprate sampling rate.
!      1.4) Estimating the unitary traction vector at each 
!           staion and component.
!   
!----------------------------------------------------------------------
!    Hugo S. Sanchez Reyes 28/11/17
!    Universite de Grenoble Alpes "UGA"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, L. Metivier and J. Virieux
!---------------------------------------------------------------------


      program main_preprocess

!----------------------Definition of variables--------------------------------------
      !COMMON VARIABLES
      implicit none
      !Define all the variables needed to read stresses
      ! and calculate the tractions associated "proc's
      ! functions", variables in include/proc.h
      INCLUDE 'proc.h'
      TYPE (mesh) :: proc_mesh

!------------------------------------------------------

!     Set working directories
     
!     Directory containing input data files P_C* TUA_C* SIGMA_**_C*
      proc_mesh%dat='dat/'
!     Directory containing output data files
      proc_mesh%out='out/'

!-----------------------------------------------------

      !Prepare unitary traction vectors as needed by 
      !the inversion code INV3DKIN
 
      !Read general information and prepare stress tensor files
      call preprocess(proc_mesh)


      end program main_preprocess
