
        subroutine init_laplacian(green_mesh)

        !COMMON VARIABLES
        IMPLICIT NONE
        INCLUDE 'green.h'
        TYPE (mesh) :: green_mesh

        integer i, iunit
        real, dimension(:), allocatable :: vector_in

        allocate(vector_in(green_mesh%msub))

        green_mesh%la(:,:) = 0.

        iunit = 32

        open(iunit,file='dat/lap.dat',status='old',action='read')
        !Five point laplacian operator [0 1 0;1 -4 1;0 1 0]
        do i=1,green_mesh%msub
         read(iunit,*) vector_in(1:green_mesh%msub)
         green_mesh%la(i,:) = vector_in(:)
        enddo

        close(iunit)

        deallocate(vector_in)
        endsubroutine init_laplacian
