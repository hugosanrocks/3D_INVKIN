      subroutine read_grad(green_mesh)

      implicit none
      ! Define all the variables needed to read models and
      ! associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      ! Variables needed only here
      real, dimension(:), allocatable :: g
      integer :: n, i, j, k, l, m

      ! ATTENTION: IF YOU WANT TO FILTER THE GRADIENT   !
      ! THAT FILTERING MUST BE DONE HERE AND NOT AFTER  !

      ! Gradient time interpolation
      call grad_time_interp(green_mesh)  

      if (green_mesh%rake_opt .eq. 1) then

       allocate(g(green_mesh%interp_i*green_mesh%msub))
       ! Flush array
       g(:)=0.


       l = 1
       j = 1
       do i=1,green_mesh%msub
        !Decompose slip vector along stk and along dip = 2 directions
        l = (i-1)*green_mesh%ncomp + 1
        do k=1,green_mesh%interp_i
         g(j)= ( green_mesh%tracint(k,l)*green_mesh%vslip(green_mesh%vnorm_i(i),1) &
  &       + green_mesh%tracint(k,l+1)*green_mesh%vslip(green_mesh%vnorm_i(i),2) &
  &       + green_mesh%tracint(k,l+2)*green_mesh%vslip(green_mesh%vnorm_i(i),3) ) 
          j = j+1
        enddo
       enddo

         green_mesh%grad1(:) = g(:)

         deallocate(g)

         !Filter gradient
         if (green_mesh%filtgrad .eq. 1) then
          call arrange_grad1d(green_mesh)
          !call smooth_lap(green_mesh)
         endif

      elseif (green_mesh%rake_opt .eq. 2) then

        allocate(g(green_mesh%interp_i*green_mesh%ncomp*green_mesh%msub))
        ! Flush array
        g(:)=0.

        !Rearrange gradient info into a 1D vector
        k = 1
        do i=1,green_mesh%msub          !number of subfaults  
          n = (i-1)*green_mesh%ncomp
          do j=1,green_mesh%ncomp       !number of components
           do m=1,green_mesh%interp_i   !number of time samples
            g(k) = green_mesh%tracint(m,j+n)
            k = k + 1
           enddo
          enddo
        enddo

        call model_d(g,green_mesh%grad2,green_mesh%interp_i,&
  &                  green_mesh%msub,green_mesh%slipm,&
  &                  green_mesh%dir_n,green_mesh%vnorm_i)

        deallocate(g)
 
        !Filter gradient
        if (green_mesh%filtgrad .eq. 1) then
         call arrange_grad(green_mesh)
        endif


      endif



      end subroutine read_grad
