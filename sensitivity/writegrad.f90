      subroutine write_grad(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      ! vecl = dimension of the problem
!     Variables needed only here
      integer :: iunit2, m, i, j, k

       !Open units where to write the solution
       iunit2 = 44
       if (green_mesh%rake_opt .eq. 1) then
          OPEN(iunit2,FILE=green_mesh%out//'kernel1d.out',&
  &       status='unknown')
       elseif (green_mesh%rake_opt .eq. 2) then
          OPEN(iunit2,FILE=green_mesh%out//'kernel2d.out',&
  &       status='unknown')
       endif

       if (green_mesh%rake_opt .eq. 1) then
         !Rearrange model info into a 1D vector
         k = 1
         do i=1,green_mesh%msub     !number of subfaults
          do m=1,green_mesh%interp_i
             write(iunit2,*) green_mesh%grad1(k)
             k = k + 1
          enddo
         enddo
       elseif (green_mesh%rake_opt .eq. 2) then
         !Rearrange model info into a 1D vector
         k = 1
         do i=1,green_mesh%msub     !number of subfaults
           do j=1,2
            do m=1,green_mesh%interp_i
             write(iunit2,*) green_mesh%grad(k)
             k = k + 1
            enddo
           enddo
         enddo
       endif

      close(iunit2)
      end subroutine write_grad
