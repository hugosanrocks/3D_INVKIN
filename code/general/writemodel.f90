      subroutine write_model(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      ! vecl = dimension of the problem
!     Variables needed only here
      integer :: iunit2, cont, cont2,  m, i, j, k, l

      !Prepare array with 3 components of slip-rate (x,y,z)
      call build_model(green_mesh)

       !Open units where to write the solution
       iunit2 = 44
       if (green_mesh%rake_opt .eq. 1) then
          OPEN(iunit2,FILE=green_mesh%dat//'modelpri1d.dat',&
  &       status='unknown')
       elseif (green_mesh%rake_opt .eq. 2) then
          OPEN(iunit2,FILE=green_mesh%dat//'modelpri.dat',&
  &       status='unknown')
       endif

       green_mesh%slipr(:,:) = 0.

       if (green_mesh%rake_opt .eq. 1) then
         !Rearrange model info into a 1D vector
         l = 1
         do i=1,green_mesh%msub     !number of subfaults
           do m=1,green_mesh%interp_i
             write(iunit2,*) green_mesh%model1(l)
             l = l + 1
           enddo
         enddo
       elseif (green_mesh%rake_opt .eq. 2) then
         !Rearrange model info into a 1D vector
         k = 1
         do i=1,green_mesh%msub     !number of subfaults
           do m=1,green_mesh%interp_i
            cont = (i-1)*(green_mesh%interp_i*2)+ m
            cont2 = cont+green_mesh%interp_i
            write(iunit2,*) green_mesh%model2(cont), &
  &                         green_mesh%model2(cont2)
           enddo
         enddo
       endif

       close(iunit2)


      end subroutine write_model



      subroutine read_modelpri(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      integer :: iunit, iunit2, m, i, j, k, cont, cont2
      integer :: samples, jump

      !call coor_trans2(green_mesh)
      green_mesh%model(:) = 0.
      green_mesh%model2(:) = 0.

      !for the first stage of the PIS (first time window)
      if (green_mesh%synwin .eq. 2) then
       samples = green_mesh%samp4win(green_mesh%synwin)
      !for the rest of the time windows
      else
       samples = green_mesh%samp4win(green_mesh%synwin-1)+green_mesh%mext
      endif

       !read solution from previous stages
       iunit=22
       OPEN(iunit,FILE=green_mesh%dat//'modelpri.dat',&
  &         status='old',action='read')

       cont = 0
       !Rearrange model info into a 1D vector
       k = 1
       do i=1,green_mesh%msub     !number of subfaults
        do m=1,samples
          cont = (i-1)*(green_mesh%interp_i*2)+ m
          cont2 = cont+green_mesh%interp_i
          read(iunit,*) green_mesh%model2(cont), &
                         green_mesh%model2(cont2)
        enddo
       enddo
       close(iunit)

       call build_model(green_mesh)


      end subroutine read_modelpri




      subroutine read_modelpri1d(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh
      integer :: iunit, i, j, k, l
      real,dimension(:),allocatable :: model1d
      integer :: samples, jump

      !call coor_trans3(green_mesh)
      green_mesh%model(:) = 0.
      green_mesh%model1(:) = 0.
      green_mesh%slipmod(:,:) = 0.

      !for the first stage of the PIS (first time window)
      if (green_mesh%synwin .eq. 2) then
       samples = green_mesh%samp4win(green_mesh%synwin)
      !for the rest of the time windows
      else
       samples = green_mesh%samp4win(green_mesh%synwin-1)+green_mesh%mext
      endif

      jump = green_mesh%interp_i - samples

      allocate(model1d(samples*green_mesh%msub))

      !Read solution from previous time windows
       iunit=22
       OPEN(iunit,FILE=green_mesh%dat//'modelpri1d.dat',&
  &         status='old',action='read')
       read(iunit,*) model1d(:)
       close(iunit)

       !Rearrange model info into a 1D vector
       j = 1
       l = 1
       do i=1,green_mesh%msub
          !Decompose slip vector along stk and along dip = 2 directions
         do k=1,samples
            green_mesh%model1(l) = model1d(j)
            j = j + 1
            l = l + 1
         enddo
         !Fill with zeros to increase model size !
         if (green_mesh%optm .eq. 2) then
            l = l + jump
         endif
         !===== Prior information option ========!
       enddo

!        if ((green_mesh%synwin .ge. 5) .and. (green_mesh%synwin .lt. 7) ) then
      !   call model_norm(green_mesh)
!print *, 'antes'
      !   call detect_motion(green_mesh)
      !   call model_norm(green_mesh)
!print *, 'despues'

!        endif


       call build_model(green_mesh)

      deallocate(model1d)
      end subroutine read_modelpri1d



      subroutine build_model(green_mesh)

      !===================================================!
      !This subroutine builds 1D array to store the three !
      !components of the slip-rate in a way that the      !
      !forward problem can be computed                    !
      !=======+===========================================!

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      integer :: i, j, k, l
      real :: vec3(1,3)

      green_mesh%model(:) = 0.

      !Choose if the model comes from  a 1D or 2D inversion
      !1D = rake_opt .eq. 1         2D = rake_opt .eq. 2
       if (green_mesh%rake_opt .eq. 1) then
         !Project to slip direction
         l = 1
         j = 1
         do i=1,green_mesh%msub
          !Decompose slip vector along stk and along dip = 2 directions
          do k=1,green_mesh%interp_i
            vec3(1,:) = green_mesh%model1(j) * green_mesh%vslip(green_mesh%vnorm_i(i),:)
            green_mesh%model(l) = vec3(1,1)
            green_mesh%model(l+green_mesh%interp_i) = vec3(1,2)
            green_mesh%model(l+green_mesh%interp_i*2) = vec3(1,3)
            l = l + 1
            j = j + 1
          enddo
          l = l+green_mesh%interp_i*2
         enddo
      elseif (green_mesh%rake_opt .eq. 2) then
        call model_c(green_mesh%model2,green_mesh%model,&
  &          green_mesh%interp_i,green_mesh%msub,green_mesh%slipm, &
  &          green_mesh%dir_n,green_mesh%vnorm_i)
      endif


      end subroutine build_model




      subroutine write_prior(green_mesh)

      !===================================================!
      !Subroutine to write prior model at current stage   !
      !of the inversion                                   !
      !===================================================!

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      integer :: i, j, k, iunit

      iunit=22

      open(iunit,file=green_mesh%dat//'current_prior.dat',&
  &        status='unknown',action='write')

      !If rake angle is assumed
      if ( green_mesh%rake_opt .eq. 1 ) then
       do i=1,green_mesh%msub
        do j=1,green_mesh%interp_i
         write(iunit,*) green_mesh%p_model1d(j,i)
        enddo
       enddo
      !if rake angle is not assumed
      elseif ( green_mesh%rake_opt .eq. 2) then
       do i=1,green_mesh%msub
        do j=1,green_mesh%interp_i
         write(iunit,*) green_mesh%p_model2d(j,i),&
  &       green_mesh%p_model2d(j,i+green_mesh%msub)
        enddo
       enddo
      endif
      close(iunit)

      endsubroutine write_prior
