      subroutine rakeout(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

      ! Arrays to save rake angle and rake angle prior
      real, dimension(:), allocatable :: angle, series
      integer, dimension(:), allocatable :: ind
      real :: val
      integer i, j, k, cont
      integer :: iunit, ii, cont2
      real pi2, trigg
      real pi/3.1415926535897932/
      pi2 = pi * 2.

      trigg = 0.1

        iunit=32
        open(iunit,file='angle.out',status='old',action='read')
 
        !samples, subfaults
        !Arrays to store the rake angles
        allocate(angle(green_mesh%interp_i*green_mesh%msub))
        allocate(ind(green_mesh%interp_i*green_mesh%msub))

        angle(:) = 0.

        cont2 = 0
        j = 1
        do i=1,green_mesh%msub
          do k=1,green_mesh%interp_i
            read(iunit,*) ii, cont, angle(j)
             if (cont .eq. 1) then
              cont2 = cont2 + 1
              ind(cont2) = ii
             endif
             j = j + 1
            enddo
            allocate(series(cont2))
            do k = 1, cont2
             series(k) = angle(ind(k))
             if (series(k) .gt. 300.) then
              series(k) = series(k) - 360.
             endif
             print *, series(k), ind(k)
            enddo
 stop
            deallocate(series)
        enddo        

        deallocate(ind)
        deallocate(angle)
      endsubroutine rakeout



