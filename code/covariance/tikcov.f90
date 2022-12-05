         program model

         !IMPLICIT NONE
         !INCLUDE 'green.h'
         !TYPE(mesh) :: green_mesh

         implicit none
         integer m, n, i, j, mm
         real vec(350,1), res1(1,350), res(1,1), res2(350,1)
         real al, be, mat(350,350), vec2(1,350)

         al = 1.d0
         be = 0.d0
         mm = 350

         do i=1,mm
          do j=1,mm
           mat(i,j) = real(j)
          enddo
          !print *, mat(i,:)
          vec(i,1) = i + (i*2)
         enddo
         !print *, 'vec', vec
         !Flush the array
         !res(:,:) = 0.d0
         !vec2(1,:) = vec(:,1)

            !Multiplication of matrices 
            call sgemm('N','N',mm,1,mm,al,        &
      &       mat,mm,vec,mm,be,res2,mm)

            !Multiplication of matrices 
            call sgemm('T','N',1,mm,mm,al,        &
      &       vec,mm,mat,mm,be,res1,1)

            !Multiplication of matrices 
            call sgemm('N','N',1,1,mm,al,        &
      &       res1,1,vec,mm,be,res,mm)
            print *, 'good', res
         
            !Multiplication of matrices 
            call sgemm('T','N',1,mm,mm,al,        &
      &       vec,1,res2,mm,be,res,1)
            print *, 'trial', res
 

         !print *, 'res2', res2
         !print *, 'res1', res1
         !print *, 'res', res
         !used only to check
         !do i=1,green_mesh%msub
         !write(33,*) green_mesh%slipr2(50,i)
         !enddo




         endprogram model

