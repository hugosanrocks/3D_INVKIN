      subroutine new_grad(green_mesh,optim)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      include 'optim_type.h'
      type (optim_type) :: optim



!     Variables needed only here
      integer :: iunit         ! dimension of the problem

       call build_model(green_mesh)

       call forward(green_mesh)

       !Compute residuals for all stations and components
       call residual(green_mesh)

         iunit=38
         open(iunit,file=green_mesh%out//'data_norm.out',status='old',access='append')
         write(iunit,*) green_mesh%costa!, optim%cpt_iter
         close(iunit)

        !TIME DOMAIN GRADIENT COMPUTATION
        call conadjtime(green_mesh)

        !build gradient in 2D and filter it
        call read_grad(green_mesh)

       !ADDITIONAL MODEL TERMS
       call penalty_terms(green_mesh)

      ! Write used to check
      ! do i=1,green_mesh%modelsize2
      !   write(81,*) green_mesh%model2(i)
      !   write(83,*) green_mesh%gradad(i)
      ! enddo

      end subroutine new_grad





      subroutine penalty_terms(green_mesh)

      !===================================================!
      !This subroutine computes the penalty terms added to!
      !the gradient either for 1D or 2D inversion. It also!
      !adds to the gradient the information required to   !
      !penalize.                                          !
      !=======+===========================================!

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer i, iunit


       !Choose if the inversion is 1D or 2D
       !1D = rake_opt .eq. 1         2D = rake_opt .eq. 2
       if (green_mesh%rake_opt .eq. 1) then

!          !Rupture time regularization 1D
!          call modeltimer1d(green_mesh)
!          print *, 'Cost: ', green_mesh%costa, 'Cost time: ',green_mesh%costm
!          green_mesh%costa = green_mesh%costa + &
!  &       green_mesh%lam1*green_mesh%costm
!          green_mesh%grad1(:) = green_mesh%grad1(:) + &
!  &       green_mesh%lam1*green_mesh%gradad1(:)

          if (green_mesh%quota1 .gt. 0.) then
           !Rupture time regularization term
           call modeltimer1d(green_mesh)
           print *, ' Misfit:', green_mesh%costa, ' Rup time misfit:',green_mesh%costm
           !Uncomment next lines to use variable regularization
           !green_mesh%balance1 = green_mesh%costm / green_mesh%costa
           !if (green_mesh%balance1 .eq. 0.) then
           ! green_mesh%lam1 = 1.
           !elseif (green_mesh%costm .gt. 1.) then
           ! green_mesh%lam1 = (1./green_mesh%balance1)*green_mesh%quota1
           !endif
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota1*green_mesh%costm
           green_mesh%grad1(:) = green_mesh%grad1(:) + &
  &        green_mesh%quota1*green_mesh%gradad1(:)
           print *, ' Cummulative misfit', green_mesh%costa
          endif


!          !Slip edge regularization 1D
!          call model_edge1d(green_mesh)
!          print *, 'Cost: ', green_mesh%costa, 'Cost edge: ',green_mesh%costm
!          green_mesh%costa = green_mesh%costa + &
!  &       green_mesh%lam2*green_mesh%costm
!          green_mesh%grad1(:) = green_mesh%grad1(:) + &
!  &       green_mesh%lam2*green_mesh%gradad1(:)

          if (green_mesh%quota2 .gt. 0.) then
           !Minimum slip at edges term
           call model_edge1d(green_mesh)  !edge effect
           print *, ' Misfit:', green_mesh%costa, ' Edge misfit:',green_mesh%costm
           !green_mesh%balance2 = green_mesh%costm / green_mesh%costa
           !if (green_mesh%balance2 .eq. 0.) then
           ! green_mesh%lam2 = 1.
           !elseif (green_mesh%costm .gt. 1.) then
           ! green_mesh%lam2 = (1./green_mesh%balance2)*green_mesh%quota2
           !endif
           !print *, green_mesh%lam2*green_mesh%balance2, 'ratio2'
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota2*green_mesh%costm
           green_mesh%grad1(:) = green_mesh%grad1(:) + &
  &        green_mesh%quota2*green_mesh%gradad1(:)
           print *, ' Cummulative misfit:', green_mesh%costa
          endif

          !Tikhonov term (model smoothing)
          if (green_mesh%quota3 .gt. 0.) then
           call model_tiko1d(green_mesh)
           print *, ' Misfit:', green_mesh%costa, ' Tikho misfit:',green_mesh%costm
           !Place to implement evolutive smoothing!
           !--------------------------------------!
           !Add misfit from model term to cost and gradient
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota3*green_mesh%costm
           green_mesh%grad1(:) = green_mesh%grad1(:) + &
  &        green_mesh%quota3*green_mesh%gradad1(:)
           print *, ' Cummulative misfit:', green_mesh%costa
          endif

          if (green_mesh%quota4 .gt. 0.) then
           !Prior model term
           call model_pri1d(green_mesh)
           !Balance the contribution to the misfit
!           green_mesh%balance4 = green_mesh%costm / green_mesh%costa
          ! if (green_mesh%costm .lt. 0.8) then
          !   green_mesh%balance4 = 1.
          ! elseif (green_mesh%costm .ge. 0.8) then
          !   green_mesh%quota4 = 0.2
          !   green_mesh%balance4 = green_mesh%costm / green_mesh%costa
          ! endif
          ! green_mesh%lam4 = (1./green_mesh%balance4)*green_mesh%quota4
           print *, ' epsilon for prior term', green_mesh%lam4

           !Add misfit from model term to cost and gradient
           print *, ' Misfit:', green_mesh%costa, ' Prior model misfit:',green_mesh%costm
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%lam4*green_mesh%costm
           green_mesh%grad1(:) = green_mesh%grad1(:) + &
  &        green_mesh%lam4*green_mesh%gradad1(:)
           !do i=1,green_mesh%modelsize1
           ! write(55,*) green_mesh%gradad(i)
           !enddo
          endif
          print *, ' Total misfit:', green_mesh%costa


       elseif (green_mesh%rake_opt .eq. 2) then
!print *, green_mesh%grad2(9295), 'inside if'
          if (green_mesh%quota1 .gt. 0.) then
           !Rupture time regularization term
           call modeltimer(green_mesh)
           print *, ' Misfit:', green_mesh%costa, ' Rup time misfit: ',green_mesh%costm
           !Uncomment next lines to use variable regularization
           !green_mesh%balance1 = green_mesh%costm / green_mesh%costa
           !if (green_mesh%balance1 .eq. 0.) then
           ! green_mesh%lam1 = 1.
           !elseif (green_mesh%costm .gt. 1.) then
           ! green_mesh%lam1 = (1./green_mesh%balance1)*green_mesh%quota1
           !endif
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota1*green_mesh%costm
           green_mesh%grad2(:) = green_mesh%grad2(:) + &
  &        green_mesh%quota1*green_mesh%gradad(:)
           print *, ' Cummulative misfit:', green_mesh%costa
          endif

          if (green_mesh%quota2 .gt. 0.) then
           !Minimum slip at edges term
           call model_edge(green_mesh)  !edge effect
           print *, ' Misfit: ', green_mesh%costa, ' Edge misfit: ',green_mesh%costm
           !green_mesh%balance2 = green_mesh%costm / green_mesh%costa
           !if (green_mesh%balance2 .eq. 0.) then
           ! green_mesh%lam2 = 1.
           !elseif (green_mesh%costm .gt. 1.) then
           ! green_mesh%lam2 = (1./green_mesh%balance2)*green_mesh%quota2
           !endif
           !print *, green_mesh%lam2*green_mesh%balance2, 'ratio2'
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota2*green_mesh%costm
           green_mesh%grad2(:) = green_mesh%grad2(:) + &
  &        green_mesh%quota2*green_mesh%gradad(:)
           print *, ' Cummulative misfit:', green_mesh%costa
          endif

          !Tikhonov term (model smoothing)
          if (green_mesh%quota3 .gt. 0.) then
           call model_tiko(green_mesh)
           print *, ' Misfit:', green_mesh%costa, ' Tikho misfit: ',green_mesh%costm
           !Place to implement evolutive smoothing!
           !--------------------------------------!
           !Add misfit from model term to cost and gradient
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%quota3*green_mesh%costm
           green_mesh%grad2(:) = green_mesh%grad2(:) + &
  &        green_mesh%quota3*green_mesh%gradad(:)
           print *, ' Cummulative misfit:', green_mesh%costa
          endif

          if (green_mesh%quota4 .gt. 0.) then
           !Prior model term
           call model_pri(green_mesh)
           !Balance the contribution to the misfit
           !if (green_mesh%costm .lt. 2.) then
           !  green_mesh%balance4 = 1.
           !print *, 'balance4'
           !else
           !  green_mesh%balance4 = green_mesh%costm / green_mesh%costa
           !endif
           !green_mesh%lam4 = (1./green_mesh%balance4)*green_mesh%quota4

           print *, ' epsilon for prior term', green_mesh%lam4

           !Add misfit from model term to cost and gradient
           print *, ' Misfit: ', green_mesh%costa, ' Prior model misfit:',green_mesh%costm
           green_mesh%costa = green_mesh%costa + &
  &        green_mesh%lam4*green_mesh%costm
           green_mesh%grad2(:) = green_mesh%grad2(:) + &
  &        green_mesh%lam4*green_mesh%gradad(:)
           !do i=1,green_mesh%modelsize2
           ! write(55,*) green_mesh%gradad(i)
           !enddo
          endif
          print *, ' Total misfit:', green_mesh%costa


!          call rake_angle(green_mesh)
         ! call rake_uni(green_mesh)
         ! print *, 'cost ang:', green_mesh%costm
          !print *, green_mesh%grad2(10), 'grad'
!          green_mesh%grad2(:) = green_mesh%grad2(:) + &
! &        green_mesh%lam5*green_mesh%gradad(:)
          !print *, green_mesh%gradad(10), 'rake grad'
          !print *, green_mesh%grad2(10), 'grad + rake'
!          green_mesh%costa = green_mesh%costa + &
!  &       green_mesh%lam5*green_mesh%costm
         ! print *, 'final cost', green_mesh%costa

          !print *, green_mesh%grad2(9295), 'af stuff rare'

       endif


      end subroutine penalty_terms
