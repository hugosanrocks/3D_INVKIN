       subroutine rotate_z(x)

       !Rotation around Z axis

       implicit none
       !Vector to rotate 90 degrees clockwise
       real, intent(inout) :: x(3)
       !Rotation matrix and dummy vector
       real conma(3,3), vec(3)
       real pi/3.1415926535897932/
       !Counters
       integer i, j

         !Rotation matrix around Z from X->N, Y->E to X->E, Y->S
         conma(1,1)=cos(-1.*pi/2.)
         conma(1,2)=-sin(-1.*pi/2.)
         conma(1,3)=0.
         conma(2,1)=sin(-1.*pi/2.)
         conma(2,2)=cos(-1.*pi/2.)
         conma(2,3)=0.
         conma(3,1)=0.
         conma(3,2)=0.
         conma(3,3)=1.

         vec(:)=0.
         do i=1,3
          do j=1,3
           vec(i)=vec(i) +conma(i,j)*x(j)
          enddo
         enddo

         !Re-asign the rotated vector
         x(:) = vec(:)

       end subroutine rotate_z
