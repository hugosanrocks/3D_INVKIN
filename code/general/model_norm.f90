      subroutine model_norm(green_mesh,cont)

      !Use this subroutine only if you know the 
      !solution you want to reconstruct.


        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
         INTEGER i, j, k, l, iunit, cont
         REAL pi/3.1415926535897932/
         REAL errorl2, errorl1

         !Initialize cumulative variable
         errorl2 = 0.
         errorl1 = 0.

         iunit=39
         open(iunit,file=green_mesh%out//'model_norm.out',status='old',access='append')

         if (green_mesh%rake_opt .eq. 1) then
          k = 1
          do j = 1,green_mesh%msub
           do i = 1,green_mesh%interp_i
            errorl2 = errorl2 + (green_mesh%target(k)- &
  &                           green_mesh%model1(k))**2.
            errorl1 = errorl1 + abs(green_mesh%target(k)- &
  &                           green_mesh%model1(k))
            !print *, error, k, green_mesh%model1(k), green_mesh%target(k)
            k = k + 1
           enddo
          enddo
         elseif (green_mesh%rake_opt .eq. 2) then
          k = 1
          do j = 1,green_mesh%msub
           do i = 1,green_mesh%interp_i*2
            errorl2 = errorl2 + (green_mesh%target(k)- &
  &                          green_mesh%model2(k))**2.
            errorl1 = errorl1 + abs(green_mesh%target(k)- &
  &                          green_mesh%model2(k))
            k = k + 1
           enddo
          enddo
         endif


         errorl2 = errorl2 * 0.5

         write(iunit,*) errorl2, cont, errorl1
         close(iunit)

         print *, errorl2, errorl1, 'model norm l2 l1'
      endsubroutine model_norm


      subroutine read_model_target(green_mesh)

      !Use this subroutine only if you know the 
      !solution you want to reconstruct.

        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
         INTEGER i, j, k, l, m, iunit
         REAL,DIMENSION(:),ALLOCATABLE :: array
         REAL pi/3.1415926535897932/
         
         allocate(array(green_mesh%msub))
         !Flush the array
         green_mesh%target(:) = 0.

         iunit=40
         open(iunit,file='dat/model_target.dat',status='old')

         if (green_mesh%rake_opt .eq. 1) then
           do i = 1,green_mesh%interp_i
             read(iunit,*) array(:)
             do k = 1,green_mesh%msub
              l = i + (k-1)* green_mesh%interp_i
              green_mesh%target(l) = array(k)
             enddo 
          enddo
         elseif (green_mesh%rake_opt .eq. 2) then
           do i = 1,green_mesh%interp_i
             read(iunit,*) array(:)
             do k=1,green_mesh%msub
              l = i + (k-1)* green_mesh%interp_i*2
              m = l + green_mesh%interp_i
              green_mesh%target(l) = array(k)* & 
                     green_mesh%vslip2(green_mesh%vnorm_i(k),1)
              green_mesh%target(m) = array(k)* &
                     green_mesh%vslip2(green_mesh%vnorm_i(k),2)
            enddo
           enddo
         endif

         

         close(iunit)

      deallocate(array)
      endsubroutine read_model_target

