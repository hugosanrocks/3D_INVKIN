         subroutine model_pri(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem(1,2), i, j, k, l, n1, n2, p
         real vec2(1,2)
         real, dimension(:,:), allocatable :: residual

         allocate(residual(green_mesh%interp_i,green_mesh%msub*2))

         mem(1,1) = 1
         mem(1,2) = 1+green_mesh%interp_i

         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            vec2(1,1) = green_mesh%model2(n1)
            vec2(1,2) = green_mesh%model2(n2)
            green_mesh%slipr2(k,i) = vec2(1,1)
            green_mesh%slipr2(k,i+green_mesh%msub) = vec2(1,2)
          enddo
         enddo

         !Current model minus prior model
         do i=1,green_mesh%msub*2
          residual(:,i) = green_mesh%slipr2(:,i) - green_mesh%p_model2d(:,i)
          green_mesh%slipr2(:,i) = residual(:,i)
         enddo

         call cov_pri(green_mesh)

         !Return additional gradient to 1D array (interp_i*msub*2)
         !=======================================!
         green_mesh%gradad(:) = 0.d0
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%Interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            vec2(1,1) = green_mesh%slipr2(k,i)
            vec2(1,2) = green_mesh%slipr2(k,i+green_mesh%msub)
            green_mesh%gradad(n1) = vec2(1,1)
            green_mesh%gradad(n2) = vec2(1,2)
          enddo
         enddo

         deallocate(residual)
         endsubroutine model_pri



         subroutine model_pri1d(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem(1,2), i, j, k, l, n1, n2, p
         real vec2(1,2)
         real, dimension(:,:), allocatable :: residual

         allocate(residual(green_mesh%interp_i,green_mesh%msub))

         mem(1,1) = 1

         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%interp_i*(i-1)+(k-1)
            vec2(1,1) = green_mesh%model1(n1)
            green_mesh%slipr2(k,i) = vec2(1,1)
          enddo
         enddo

         !Current model minus prior model
         do i=1,green_mesh%msub
          residual(:,i) = green_mesh%slipr2(:,i) - green_mesh%p_model1d(:,i)
          green_mesh%slipr2(:,i) = residual(:,i)
         enddo

         call cov_pri1d(green_mesh)

         !Return additional gradient to 1D array (interp_i*msub*2)
         !=======================================!
         green_mesh%gradad(:) = 0.d0
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%Interp_i*(i-1)+(k-1)
            vec2(1,1) = green_mesh%slipr2(k,i)
            green_mesh%gradad1(n1) = vec2(1,1)
          enddo
         enddo

         deallocate(residual)
         endsubroutine model_pri1d
