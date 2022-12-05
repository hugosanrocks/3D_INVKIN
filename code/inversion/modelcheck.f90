         subroutine model_check(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem(1,2), i, j, k, l, n1, n2
         real vec2(1,2)
         real matt2(3,2)

         !Flush the array
         green_mesh%slipr2(:,:) = 0.d0

         mem(1,1) = 1
         mem(1,2) = 1+green_mesh%interp_i

         !Arrange into an array (interp_i x msub*2)
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%Interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            vec2(1,1) = green_mesh%model2(n1)
            vec2(1,2) = green_mesh%model2(n2)
            green_mesh%slipr2(k,i) = vec2(1,1)
            green_mesh%slipr2(k,i+green_mesh%msub) = vec2(1,2)
          enddo
         enddo

         !do i=1,875
         ! write(32,*) green_mesh%slipr2(i,1), green_mesh%slipr2(i,162)
         !enddo

         !Fill with zeros where energy should not be
         do i=1,2
          do j=1,green_mesh%msub
           green_mesh%slipr2(1:green_mesh%rsamp(j),j+(i-1)*green_mesh%msub) = 0.
          enddo
         enddo
         !do i=1,875
         ! write(33,*) green_mesh%slipr2(i,1), green_mesh%slipr2(i,162), green_mesh%slipr2(i,289), green_mesh%slipr2(i,162+288)
         !enddo

         !Return additional model to 1D array (interp_i*msub*2)
         !=======================================!
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
           n1 = mem(1,1) + green_mesh%Interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            vec2(1,1) = green_mesh%slipr2(k,i)
            vec2(1,2) = green_mesh%slipr2(k,i+green_mesh%msub)
            green_mesh%model2(n1) = vec2(1,1)
            green_mesh%model2(n2) = vec2(1,2)
          enddo
         enddo



         endsubroutine model_check
