       ! program main

       ! real val(1,1), vect(1,3), vec2(1,2), vec3(1,3), matt(2,3), matt2(3,2), vec32(3,1)
       ! integer i, j

       ! val(1,1) = 5
       ! vec2(1,:) = [  -0.98481,   0.17365 ]
       ! vec3(1,:) = [  -0.040009,  0.992945,   0.111619]
       ! vect(1,:) = vec3(1,:)
       ! vec32(:,1) = vec3(1,:)
       ! matt(1,:) = [   0.17365,  -0.98481,   0.00000]
       ! matt(2,:) = [   0.75441,   0.13302,   0.64279]

       ! matt2(1,:) = matt(:,1)
       ! matt2(2,:) = matt(:,2)
       ! matt2(3,:) = matt(:,3)
       ! vec3(1,:) = val(1,1) * vec3(1,:)
       ! vec2(:,:) = 0.

       !  call matrix1_3(val(1,1),vec3,vect)
       !  call mat3_2(vec3,vec2,matt2)
       !  call matrix3_2(vec3,vec2,matt2)
       !  call mat2_3(vec2,vec3,matt)
       !  call matrix2_3(vec2,vec3,matt)
       !  call matrix3_1(vec3,val(1,1),vect)

       ! endprogram main


         subroutine matrix1_3(val,vec3,vect)

         real, intent(inout) :: val(1,1), vec3(1,3), vect(1,3)
         real vs3(1,1), vm(1,3), res(1,3), al, be

         vec3(:,:) = 0.
         vec3(1,:) = val(1,1) * vect(1,:)
         !print *, vec3,'fin vec'

         return
         endsubroutine matrix1_3


       subroutine mat3_2(vec3,vec2,matt2)
         real, intent(inout) :: vec3(1,3), vec2(1,2), matt2(3,2)
         integer i, j
         real al, be

         al = 1.d0
         be = 0.d0
         vec2(:,:) = 0.d0
         !print *, 'vec3 1', vec3
         call sgemm('N','N',1,2,3,al,        &
      &     vec3,1,matt2,3,be,vec2,1)
        !print *, 'vec2 mkl', vec2

       return
       endsubroutine mat3_2



         subroutine matrix3_2(vec3,vec2,matt2)

         real, intent(inout) :: vec3(1,3), vec2(1,2), matt2(3,2)
         integer i, j

         vec2(:,:) = 0.
         do i=1,2
           do j=1,3
            vec2(1,i) = vec2(1,i) + vec3(1,j)*matt2(j,i)
           enddo
          enddo
          !print *, 'vec2', vec2

         return
         endsubroutine matrix3_2



      subroutine mat2_3(vec2,vec3,matt)
       real, intent(inout) :: vec2(1,2), vec3(1,3), matt(2,3)
       integer i, j
       real al, be

       al = 1.d0
       be = 0.d0
       vec3(:,:) = 0.d0
       !print *, 'vec2 1', vec2
       call sgemm('N','N',1,3,2,al,        &
      &     vec2,1,matt,2,be,vec3,1)
       !print *, 'vec3 mkl', vec3

       return
       endsubroutine mat2_3



         subroutine matrix2_3(vec2,vec3,matt)

         real, intent(inout) :: vec2(1,2), vec3(1,3), matt(2,3)
         integer i, j

         vec3(:,:) = 0.
         do i=1,3
           do j=1,2
            vec3(1,i) = vec3(1,i) + vec2(1,j)*matt(j,i)
           enddo
         enddo
         !print *, 'vec3', vec3

         return
         endsubroutine matrix2_3


         subroutine matrix3_1(vec3,val,vect)

         real, intent(inout) :: vec3(1,3), val(1,1), vect(1,3)
         integer i

         val(1,1) = 0.
         do i=1,3
          val(1,1) = val(1,1) + vec3(1,i)*vect(1,i)
         enddo
         !print *, 'val', val

         return
         endsubroutine matrix3_1



