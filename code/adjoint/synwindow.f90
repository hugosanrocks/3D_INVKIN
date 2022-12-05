      subroutine windows(green_mesh)

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer i, iunit, j
      real, dimension(:), allocatable :: vector_in


      allocate(vector_in(green_mesh%wininv))
      allocate(green_mesh%twin(green_mesh%nsta,green_mesh%wininv))
      call monitor_alloc('66%twin   ',' allo')


      iunit = 22
      open(iunit,file=green_mesh%dat//'windows.info',status='old',&
  &        action='read')
      do i= 1,green_mesh%nsta
        read(iunit,*) j, vector_in(1:green_mesh%wininv)
        green_mesh%twin(i,:) = vector_in(1:green_mesh%wininv)
        green_mesh%samwin(i,:) = int(green_mesh%twin(i,:)/real(green_mesh%intdt)) + 1
      enddo

      close(iunit)

      deallocate(vector_in)
      endsubroutine windows
