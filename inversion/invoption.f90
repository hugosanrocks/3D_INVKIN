     subroutine invoption(green_mesh)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      INTEGER i      
      REAL start, fini

      !===================================================================!
      ! This subroutine controls the type of inversion to be performed.   !
      ! The rake angle can be set as a known or unknown value. (rake_opt) !          !
      ! The inversion can be done for the whole rupture time or for       !
      ! progressive space-time windows                          (for_opt) !
      !===================================================================! 

      !===================================================================!
      !                   green_mesh%rake_opt                             !
      ! inversion options for assuming a known or unknown rake            !
      ! 1D = rake_opt .eq. 1         2D = rake_opt .eq. 2                 !
      !===================================================================!

      !===================================================================!
      !                   green_mesh%for_opt                              !
      ! inversion options for initial model and time windows              !
      ! 1 = start from slip-rate 0 and invert the whole rupture history   !
      ! 2 = start from 0 and invert progressively space-time windows      !
      !     using previously estimated windows                            !
      !     as initial model for the next one                             !
      ! 3 = start from a given initial model and invert the whole history !
      !===================================================================!

      !===================================================================!
      !Initialize first values of hyper-parameters !
      ! green_mesh%lam1 = 1.
      ! green_mesh%lam2 = 1.
       !green_mesh%lam3 = 1.
       !green_mesh%lam4 = 1.
      !===================================================================!
      !Real weigths are control by green_mesh%quota*
      !Defined at dat/fwioption.info
      !===================================================================!

      !Rake is fixed
      if (green_mesh%rake_opt .eq. 1) then
         write(*,*) ' Rake is not an unknown rake_opt=1'
         if (green_mesh%for_opt .eq. 1) then
           write(*,*) ' Inversion option 1 (initial model 0) '
           call initwin(green_mesh)
           !Prepare prior model
           call read_prior_model(green_mesh)
           !call prepare_preconditioner(green_mesh)
           call exp_covar(green_mesh)
           call coor_trans(green_mesh)
           print *, ' Initial model read from dat/vitesse.out'
           !if target solution is known
           call read_model_target(green_mesh)
           i=-1
           call model_norm(green_mesh,i)
         elseif (green_mesh%for_opt .eq. 2) then
           write(*,*) ' Inversion option 2 (time windows) '
           call cpu_time(start)
           do i = 1,green_mesh%wininv - 1
             !green_mesh%lam4 = green_mesh%lam4+0.000001
             green_mesh%dowin = i
             green_mesh%synwin = i + 1
             call initwin(green_mesh)
             !Prepare prior model
             call read_prior_model(green_mesh)
             call exp_covar(green_mesh)
             call read_modelpri1d(green_mesh)
             print *, ' Initial model read from dat/modelpri.dat'
             !if target solution is known
             call read_model_target(green_mesh)
             call model_norm(green_mesh,i)
             !Number of windows growing the model space
             green_mesh%flag_d = 'PRED'
             call prediction(green_mesh)

             call destroywin(green_mesh)
           enddo
           call cpu_time(fini)
           WRITE(6, *) '==============================================='
            WRITE(6,*)  ' Time used for inversion: ', fini-start,'(sec)'
           WRITE(6, *) '==============================================='
           green_mesh%synwin = green_mesh%wininv-1
           green_mesh%mext = 0
           write(*,*) ' Final iterations and focus of energy'
           call initwin(green_mesh)
           !Prepare prior model
           call read_prior_model(green_mesh)
           call exp_covar(green_mesh)
           green_mesh%flag_d = 'DETC'
           call detect(green_mesh)
           write(*,*) ' Started final model and regularization '
         !Whole rupture time inversion from dat/modelpri.out
         endif
      !Rake is not fixed
      elseif (green_mesh%rake_opt .eq. 2) then
         write(*,*) ' Rake is an unknown rake_opt=2'
         ! Whole rupture time inversion from dat/vitesse.out
         if (green_mesh%for_opt .eq. 1) then
           write(*,*) ' Inversion option 1 (initial model 0) '
           call initwin(green_mesh)
           !Prepare prior model
           call read_prior_model(green_mesh)
           call exp_covar(green_mesh)
           call coor_trans(green_mesh)
           print *, ' Initial model read from dat/vitesse.out'
           !if target solution is known
           call read_model_target(green_mesh)
           i=-1
           call model_norm(green_mesh,i) 
        !Progressive time-space inversion
           elseif (green_mesh%for_opt .eq. 2) then
           write(*,*) ' Inversion option 2 (time windows) '
           do i = 1,green_mesh%wininv - 1
             green_mesh%dowin = i
             green_mesh%synwin = i + 1
             call initwin(green_mesh)
             !Prepare prior model
             call read_prior_model(green_mesh)
             call exp_covar(green_mesh)
             call read_modelpri(green_mesh)
             print *, ' Initial model read from dat/modelpri.dat'
             !if target solution is known
             call read_model_target(green_mesh)
             call model_norm(green_mesh,i)
             !Number of windows growing the model space
             green_mesh%flag_d = 'PRED'
             call prediction(green_mesh)
             call rake_uni(green_mesh)
             call destroywin(green_mesh)
           enddo
           call cpu_time(fini)
           WRITE(6, *) '==============================================='
            WRITE(6,*)  ' Time used for inversion: ', fini-start,'(sec)'
           WRITE(6, *) '==============================================='
           green_mesh%synwin = green_mesh%wininv-1
           green_mesh%mext = 0
           write(*,*) ' Final iterations and focus of energy'
           call initwin(green_mesh)
           !Prepare prior model
           call read_prior_model(green_mesh)
           call exp_covar(green_mesh)
           green_mesh%flag_d = 'DETC'
           call detect(green_mesh)
           write(*,*) ' Started final model and regularization '
         !Whole rupture time inversion from dat/modelpri.out
         elseif (green_mesh%for_opt .eq. 3) then
           write(*,*) ' Inversion option 3 (initial model from file) '
           call initwin(green_mesh)
           call exp_covar(green_mesh)
           green_mesh%flag_d = 'DETC'
           green_mesh%synwin = 8
           green_mesh%mext = 0
           call detect(green_mesh)
           print *, ' Initail model from dat/modelpri.dat'
         else
           print *, ' Wrong option (for_opt) check in  dat/syn.info'
         endif
      else
        write(*,*) ' Wrong option (rake_opt) check in dat/focal.info'
      endif

     !if whole recordings is used
     green_mesh%synwin = green_mesh%wininv

     endsubroutine invoption
