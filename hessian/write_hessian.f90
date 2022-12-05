      subroutine write_hessian(green_mesh)


      IMPLICIT NONE
      ! Define all the GLOBAL variables needed
      ! variables inside include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only by this program

      integer i, j, k, iunit, msize


      msize = green_mesh%modelsize1

      iunit = 81
      open(iunit,file='hessian_column.bin',status='unknown',&
  &        action='write',form='unformatted',access='direct',&
  &        recl=msize)

      j = green_mesh%hess_cont
      write(iunit,rec=j) green_mesh%grad1
      close(iunit)

      !do j=1,green_mesh%modelsize1
      ! write(55,*) green_mesh%grad1(j)
      !enddo

      endsubroutine write_hessian


 



