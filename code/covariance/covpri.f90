         subroutine cov_pri(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh


         integer m, n, i, j, mm, iunit, sub, k, kk
         real cres(1,1)
         real al, be, cost
         real,dimension(:,:),allocatable :: vec, vec2, res, matt

         !require memory
         allocate(matt(green_mesh%msub,green_mesh%msub))
         allocate(vec2(1,green_mesh%msub))
         allocate(vec(green_mesh%msub,1))
         allocate(res(green_mesh%msub,1))

         !variables used by MKL library
         !for matrix multiplication
         !C = alphaA*B + betaC
         al = 1.d0
         be = 0.d0
         mm = green_mesh%msub
         green_mesh%costm = 0.
         sub=green_mesh%msub

         !Flush the array
         res(:,:) = 0.d0
         matt(:,:) = 0.

         !Arrange 1 fault snapshot in a vector
         do i=1,2
          m=1+(i-1)* mm
          n=1+(i-1)* mm + mm - 1
          kk = 1
          do j=1,green_mesh%interp_i
           do k=1,green_mesh%msub
             matt(k,k) =  green_mesh%diag(kk)
             kk = kk+1
           enddo
           vec(1:green_mesh%msub,1) = green_mesh%slipr2(j,m:n)
           vec2(1,:) = vec(:,1)
            !Multiplication of matrices 
            call sgemm('N','N',mm,1,mm,al,        &
      &       matt,mm,vec,mm,be,res,mm)
             !Multiplication of matrices 
            call sgemm('N','N',1,1,mm,al,        &
      &       vec2,1,res,mm,be,cres,mm)
            green_mesh%costm = green_mesh%costm + cres(1,1)
           green_mesh%slipr2(j,m:n) = res(1:mm,1)
          enddo
         enddo

         !multiply model misfit by 1/2
         green_mesh%costm = 0.5 * green_mesh%costm

         deallocate(matt,vec,vec2,res)
         endsubroutine cov_pri




         subroutine cov_pri1d(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh


         integer m, n, i, j, mm, iunit, sub, k, kk
         real cres(1,1)
         real al, be, cost
         real,dimension(:,:),allocatable :: vec, vec2, res, matt

         !require memory
         allocate(matt(green_mesh%msub,green_mesh%msub))
         allocate(vec2(1,green_mesh%msub))
         allocate(vec(green_mesh%msub,1))
         allocate(res(green_mesh%msub,1))

         !variables used by MKL library
         !for matrix multiplication
         !C = alphaA*B + betaC
         al = 1.d0
         be = 0.d0
         mm = green_mesh%msub
         green_mesh%costm = 0.
         sub=green_mesh%msub

         !Flush the array
         res(:,:) = 0.d0
         matt(:,:) = 0.

         !Arrange 1 fault snapshot in a vector
         do i=1,1
          m=1+(i-1)* mm
          n=1+(i-1)* mm + mm - 1
          kk = 1
          do j=1,green_mesh%interp_i
!          read(iunit,rec=j) matt         !Read matrix for time regul
           do k=1,green_mesh%msub
             matt(k,k) = green_mesh%diag(kk)
             kk = kk+1
           enddo
           vec(1:green_mesh%msub,1) = green_mesh%slipr2(j,m:n)
           vec2(1,:) = vec(:,1)
            !Multiplication of matrices
            call sgemm('N','N',mm,1,mm,al,        &
      &       matt,mm,vec,mm,be,res,mm)
             !Multiplication of matrices
            call sgemm('N','N',1,1,mm,al,        &
      &       vec2,1,res,mm,be,cres,mm)
            green_mesh%costm = green_mesh%costm + cres(1,1)
           green_mesh%slipr2(j,m:n) = res(1:mm,1)
          enddo
         enddo

         !multiply model misfit by 1/2
         green_mesh%costm = 0.5 * green_mesh%costm

         deallocate(matt,vec,vec2,res)
         endsubroutine cov_pri1d
