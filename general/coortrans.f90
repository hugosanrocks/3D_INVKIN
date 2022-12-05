      subroutine coor_trans(green_mesh)

        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
         INTEGER i, j, k, l, iunit
         REAL vec3(1,3)

         real, dimension(:), allocatable :: vector_in
         !used only to avoid I/O isuues
         allocate(vector_in(green_mesh%msub))



         !SAME AS AXITRA
         !****************************************************************!
         ! - STRIKE : strike of the fault measured clockwise from North. 
         !   DIP : dip of the fault. 
         !   RAKE : direction of slip of the hanging wall relatively to 
         !   the foot wall. It is measured counterclockwise from the 
         !   strike direction. If faulting is right-lateral (like the San 
         !   Andreas or the North Anatolian), rake=0.
         !****************************************************************!

         !All rotations have been moved to normal_vectos.f90
         !Option to estimate explicitly the Hessian
         !hess_opt=1 estimate, else = do not estimate
         if (green_mesh%hess_opt .eq. 1) then

         else
           !Read the slip-rate initial model (modulus of vector)
           iunit=30
           open(iunit,file=green_mesh%dat//'vitesse.out',status='unknown',action='read')
           do i=1,green_mesh%interp_i      !slipsam = interp_i
            read(iunit,*) vector_in(1:green_mesh%msub)
            green_mesh%slipmod(i,:) = vector_in(1:green_mesh%msub)
           enddo
           close(iunit)
         endif

      if (green_mesh%rake_opt .eq. 1) then

        if (green_mesh%hess_opt .eq. 1) then
          write(*,*) ' Hessian option 1 turned on '
          write(*,*) ' Setting cannonical vectors '
          green_mesh%model1(:) = 0.
          green_mesh%model1(green_mesh%hess_cont) = 1.
        elseif (green_mesh%hess_opt .eq. 2) then
          write(*,*) ' Hessian option 2 turned on '
          write(*,*) ' Product A^{T}b'
          green_mesh%model1(:) = 0.
        else
          write(*,*) ' Hessian option other turned on '
          write(*,*) ' No cannonical vectors '
          !Arrange slip rate model in 1D vector
          l = 1
          do i=1,green_mesh%msub
           !Order slip model into a 1D array, rake is known
           do k=1,green_mesh%interp_i
             green_mesh%model1(l) = green_mesh%slipmod(k,i)
             l = l + 1
           enddo
          enddo
        endif
        !Project to the known 2D slip direction        
        !rake is known
        l = 1
        j = 1
        do i=1,green_mesh%msub
         !Decompose slip vector along stk and along dip = 2 directions
         do k=1,green_mesh%interp_i
          vec3(1,:) = green_mesh%model1(j)*green_mesh%vslip(green_mesh%vnorm_i(i),:)
          green_mesh%model(l) = vec3(1,1)
          green_mesh%model(l+green_mesh%interp_i) = vec3(1,2)
          green_mesh%model(l+green_mesh%interp_i*2) = vec3(1,3)
          l = l + 1
          j = j + 1
         enddo
         l = l+green_mesh%interp_i*2
        enddo

      !If rake is not known, slip is a 2D vector
      elseif (green_mesh%rake_opt .eq. 2) then

        !Arrange slip rate model in 1D array, 2D directions
        l = 1
        do i=1,green_mesh%msub
        !Decompose slip vector along stk and along dip = 2 directions
         do j=1,2
          do k=1,green_mesh%interp_i
           green_mesh%model2(l) = green_mesh%slipmod(k,i)*green_mesh%vslip2(green_mesh%vnorm_i(i),j)
           l = l + 1
          enddo
         enddo
        enddo


        !Change slip vector model from along stk and dip to (x,y,z)
        call model_c(green_mesh%model2,green_mesh%model,&
  &          green_mesh%interp_i,green_mesh%msub,green_mesh%slipm,&
             green_mesh%dir_n,green_mesh%vnorm_i)



      endif

      deallocate(vector_in)
      endsubroutine coor_trans




      !subroutine prepare_prior_model(green_mesh)

        !COMMON VARIABLES
      !   IMPLICIT NONE
      !   INCLUDE 'green.h'
      !   TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
      !   INTEGER i, j, k, l, iunit, mem(1,2), n1, n2
      !   REAL val(1,1)
      !   REAL res(1,2)

         !SAME AS AXITRA
         !****************************************************************!
         ! - STRIKE : strike of the fault measured clockwise from North. 
         !   DIP : dip of the fault. 
         !   RAKE : direction of slip of the hanging wall relatively to 
         !   the foot wall. It is measured counterclockwise from the 
         !   strike direction. If faulting is right-lateral (like the San 
         !   Andreas or the North Anatolian), rake=0.
         !****************************************************************!


      !   endsubroutine prepare_prior_model


        subroutine read_prior_model(green_mesh)

        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
         INTEGER i, j, iunit
         REAL vec2d(2)

         real, dimension(:), allocatable :: vector_in

         !used only to avoid I/O isuues
         allocate(vector_in(green_mesh%msub))


         !Prior model for 1D slip-rate functions
         if (green_mesh%rake_opt .eq. 1) then
           iunit=30
           open(iunit,file=green_mesh%dat//'prior_model.dat',status='unknown',action='read')
           do i=1,green_mesh%interp_i      !slipsam = interp_i
            read(iunit,*) vector_in(1:green_mesh%msub)
            green_mesh%p_model1d(i,1:green_mesh%msub) = vector_in(1:green_mesh%msub)
           enddo
           close(iunit)
         !Prior model for 2D slip-rate functions
         elseif (green_mesh%rake_opt .eq. 2) then
           iunit=30
           open(iunit,file=green_mesh%dat//'prior_model.dat',status='unknown',action='read')
           do i=1,green_mesh%interp_i      !slipsam = interp_i
            read(iunit,*) green_mesh%p_model2d(i,1:green_mesh%msub)
            do j=1,green_mesh%msub
             !Decompose vector along stk and dip = 2 directions
             vec2d(:) = green_mesh%p_model2d(i,j)*green_mesh%vslip2(green_mesh%vnorm_i(j),:)
             green_mesh%p_model2d(i,j) = vec2d(1)
             green_mesh%p_model2d(i,j+green_mesh%msub) = vec2d(2)
            enddo
           enddo
           close(iunit)
         endif

         deallocate(vector_in)
         endsubroutine read_prior_model


      subroutine detect_focmec(green_mesh)

      implicit none
      !Define all the variables needed to read stresses and calculate the tractions
      !associated "Green's functions", variables in include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only by this program
      integer i, iunit

      iunit=25
      open(iunit,file=green_mesh%dat//'/map_focal.dat',status='old',&
  &        action='read')
      do i=1,green_mesh%msub
        read(iunit,*) green_mesh%vnorm_i(i)
      enddo
      close(iunit)

      endsubroutine detect_focmec


