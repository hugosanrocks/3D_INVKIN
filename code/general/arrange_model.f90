      subroutine arrange_model(green_mesh,model,vecl)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: i, j, k       ! counters

      integer, intent(inout) :: vecl                 ! dimension of the problem
      real, intent(inout) :: model(vecl)             ! current model & gradient vectors

          !Rearrange model and gradient info in a vector
          k=1
          do j=1,green_mesh%ncomp
           do i=1,green_mesh%interp_i
            model(k)=green_mesh%slipr(i,j)        !Arrange model in 1D vector
            k=k+1
           enddo
          enddo


         end subroutine arrange_model
