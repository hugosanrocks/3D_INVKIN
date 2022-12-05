      subroutine vectors_geometry(green_mesh)

        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE (mesh) :: green_mesh

         !Variables needed by this subroutine
         INTEGER i, j
         REAL pi/3.1415926535897932/
         !REAL v3(2), vec3(1,3), vec2(1,2), matt(2,3), matt2(3,2), vect(1,3), val(1,1)
         !REAL res(1,2), vs3(1,3), vm(3,2), al, be
         REAL vector(3)

         !SAME AS AXITRA
         !****************************************************************!
         ! - STRIKE : strike of the fault measured clockwise from North. 
         !   DIP : dip of the fault. 
         !   RAKE : direction of slip of the hanging wall relatively to 
         !   the foot wall. It is measured counterclockwise from the 
         !   strike direction. If faulting is right-lateral (like the San 
         !   Andreas or the North Anatolian), rake=0.
         !****************************************************************!

         !Add 90 degrees to rake to be consistent with AXITRA
         green_mesh%rak=green_mesh%rak+pi

         !Unitary slip vector
         !Same as AXITRA coordinate system
         green_mesh%vslip(green_mesh%dir_i,1) = (cos(green_mesh%rak)*cos(green_mesh%stk)) + &   !-1.0
  &      sin(green_mesh%rak)*cos(green_mesh%dip)*sin(green_mesh%stk)
         green_mesh%vslip(green_mesh%dir_i,2) =( cos(green_mesh%rak)*sin(green_mesh%stk)) - & 
  &      (sin(green_mesh%rak)*cos(green_mesh%dip)*cos(green_mesh%stk))
         green_mesh%vslip(green_mesh%dir_i,3) = -1.0*sin(green_mesh%rak)*sin(green_mesh%dip)  !-1 z para arriba

         !Fault unitary normal vector used for transformation (Stein p.218 eq.2)
         !Same as AXITRA coordinate system
         green_mesh%vnorm(green_mesh%dir_i,1)= -1.*sin(green_mesh%dip)*sin(green_mesh%stk)
         green_mesh%vnorm(green_mesh%dir_i,2)= sin(green_mesh%dip)*cos(green_mesh%stk)
         green_mesh%vnorm(green_mesh%dir_i,3)= -1.*cos(green_mesh%dip)


         !Unitary strike vector
         green_mesh%vstk(green_mesh%dir_i,1)=cos(green_mesh%stk)
         green_mesh%vstk(green_mesh%dir_i,2)=sin(green_mesh%stk)
         green_mesh%vstk(green_mesh%dir_i,3)=0.


      !   write(6,*) ' Slip vector    :', green_mesh%vslip(green_mesh%dir_i,:)
      !   write(6,*) ' Strike vector  :', green_mesh%vstk(green_mesh%dir_i,:)
      !   write(6,*) ' norm vector     :', green_mesh%vnorm(green_mesh%dir_i,:)



         !Rotate vectors to be consistent with mesh GEODG3D
         vector(:) = green_mesh%vnorm(green_mesh%dir_i,:)
         call rotate_z(vector)
         green_mesh%vnorm(green_mesh%dir_i,:) = vector(:)

         vector(:) = green_mesh%vslip(green_mesh%dir_i,:)
         call rotate_z(vector)
         green_mesh%vslip(green_mesh%dir_i,:) = vector(:)
         vector(:) = green_mesh%vstk(green_mesh%dir_i,:)
         call rotate_z(vector)
         green_mesh%vstk(green_mesh%dir_i,:) = vector(:)

         !Unitary dip vector
         green_mesh%vdip(green_mesh%dir_i,1)= green_mesh%vstk(green_mesh%dir_i,2)*&
  &        green_mesh%vnorm(green_mesh%dir_i,3)-green_mesh%vnorm(green_mesh%dir_i,2)*&
  &        green_mesh%vstk(green_mesh%dir_i,3)
         green_mesh%vdip(green_mesh%dir_i,2)= -1.*(green_mesh%vstk(green_mesh%dir_i,1)*&
  &        green_mesh%vnorm(green_mesh%dir_i,3)-green_mesh%vnorm(green_mesh%dir_i,1)*&
  &        green_mesh%vstk(green_mesh%dir_i,3))
         green_mesh%vdip(green_mesh%dir_i,3)= green_mesh%vstk(green_mesh%dir_i,1)*&
  &        green_mesh%vnorm(green_mesh%dir_i,2)-green_mesh%vnorm(green_mesh%dir_i,1)*&
  &        green_mesh%vstk(green_mesh%dir_i,2)

         green_mesh%slipm(green_mesh%dir_i,1,:) = green_mesh%vstk(green_mesh%dir_i,:)
         green_mesh%slipm(green_mesh%dir_i,2,:) = green_mesh%vdip(green_mesh%dir_i,:)

         !Decompose slip vector along strike and dip
         green_mesh%vslip2(green_mesh%dir_i,:)=0.
         do i=1,2
          do j=1,3
           green_mesh%vslip2(green_mesh%dir_i,i)=green_mesh%vslip2(green_mesh%dir_i,i) + &
  &             green_mesh%slipm(green_mesh%dir_i,i,j)*green_mesh%vslip(green_mesh%dir_i,j)
          enddo
         enddo

         !Print only to check decompositon along strike and dip
         vector(1:3) = green_mesh%vslip(green_mesh%dir_i,1:3)
         write(6,*) ' Slip vector    :', vector
         vector(1:3) = green_mesh%vstk(green_mesh%dir_i,1:3)
         write(6,*) ' Strike vector  :', vector
         vector(1:3) = green_mesh%vdip(green_mesh%dir_i,1:3)
         write(6,*) ' Dip vector     :', vector
         vector(1:2) = green_mesh%vslip2(green_mesh%dir_i,1:2)
         write(6,*) ' 2D Slip vector :', vector(1:2)

      endsubroutine vectors_geometry


      subroutine vector_transformation(stk,dip,rake,vector_out)

        !COMMON VARIABLES
         IMPLICIT NONE

         !Variables needed by this subroutine
         REAL, INTENT(INOUT) :: rake, stk, dip, vector_out(2)
         INTEGER i, j
         REAL pi/3.1415926535897932/
         REAL vector(3), vslip(3), vdip(3), vstk(3), vnorm(3), slipm(2,3)

         !SAME AS AXITRA
         !****************************************************************!
         ! - STRIKE : strike of the fault measured clockwise from North. 
         !   DIP : dip of the fault. 
         !   RAKE : direction of slip of the hanging wall relatively to 
         !   the foot wall. It is measured counterclockwise from the 
         !   strike direction. If faulting is right-lateral (like the San 
         !   Andreas or the North Anatolian), rake=0.
         !****************************************************************!

         !Add 90 degrees to rake to be consistent with AXITRA
         rake = rake*pi/180.
         stk = stk*pi/180.
         dip = dip*pi/180.
         rake= rake+pi

         !Unitary slip vector
         !Same as AXITRA coordinate system
         vslip(1) = (cos(rake)*cos(stk)) + &   !-1.0
  &      sin(rake)*cos(dip)*sin(stk)
         vslip(2) =( cos(rake)*sin(stk)) - & 
  &      (sin(rake)*cos(dip)*cos(stk))
         vslip(3) = -1.0*sin(rake)*sin(dip)

         !Fault unitary normal vector used for transformation (Stein p.218 eq.2)
         !Same as AXITRA coordinate system
         vnorm(1)= -1.*sin(dip)*sin(stk)
         vnorm(2)= sin(dip)*cos(stk)
         vnorm(3)= -1.*cos(dip)


         !Unitary strike vector
         vstk(1)=cos(stk)
         vstk(2)=sin(stk)
         vstk(3)=0.


      !   write(6,*) ' Slip vector    :', green_mesh%vslip(green_mesh%dir_i,:)
      !   write(6,*) ' Strike vector  :', green_mesh%vstk(green_mesh%dir_i,:)
      !   write(6,*) ' norm vector     :', green_mesh%vnorm(green_mesh%dir_i,:)



         !Rotate vectors to be consistent with mesh GEODG3D
         vector(:) = vnorm(:)
         call rotate_z(vector)
         vnorm(:) = vector(:)
         vector(:) = vslip(:)
         call rotate_z(vector)
         vslip(:) = vector(:)
         vector(:) = vstk(:)
         call rotate_z(vector)
         vstk(:) = vector(:)

         !Unitary dip vector
         vdip(1)= vstk(2)*vnorm(3)-vnorm(2)*vstk(3)
         vdip(2)= -1.*(vstk(1)*vnorm(3)-vnorm(1)*vstk(3))
         vdip(3)= vstk(1)*vnorm(2)-vnorm(1)*vstk(2)

         slipm(1,:) = vstk(:)
         slipm(2,:) = vdip(:)

         !Decompose slip vector along strike and dip
         vector_out(:)=0.
         do i=1,2
          do j=1,3
           vector_out(i)=vector_out(i) + slipm(i,j)*vslip(j)
          enddo
         enddo

         !Print only to check decompositon along strike and dip
         !write(6,*) ' Slip vector    :', vslip(:)
         !write(6,*) ' Strike vector  :', vstk(:)
         !write(6,*) ' Dip vector     :', vdip(:)
         !write(6,*) ' 2D Slip vector :', vector_out(:)

      endsubroutine vector_transformation
