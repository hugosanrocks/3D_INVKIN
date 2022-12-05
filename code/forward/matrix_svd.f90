        program matrix_svd

    implicit none

       integer ntime, ncomp, msub, stcomp
       integer lensyn, lenslip, nx, ny, j, k, iunit
       INTEGER M2, N2, i
       INTEGER lenH, lenU, lenVT
       INTEGER LWMAX, m, n, val
       INTEGER INFO, LWORK
       INTEGER,dimension(:),allocatable :: IWORK
       REAL,DIMENSION(:),ALLOCATABLE :: S, WORK, WORK2
       REAL,DIMENSION(:,:),ALLOCATABLE :: H, U, VT, SI, VSI, VSIUT, UT, V, ATb, xsol
       REAL alpha, beta
       INTEGER option

       !Used for Hessian vector product
       !Coeficient for matrix multiplication C = alpha*A*B + beta*C
       alpha = 1.
       beta = 0.

       write(*,*) 'option 1 computes     option 2 reads' 
       read(*,*) option

       !Size of matrices H = V S U, all are the same
       !size of the square matrix to decompose m2 x m2
       m2 = 189*50
       n2 = m2
       lenH = m2
       lenU = m2
       lenVT = n2
       lwmax = 1000

       !Allocate array asked by LAPACK
       allocate(IWORK(8*N2))
       !Allocate matrices to save A, U, S, VT
       allocate(H(lenH,lenH),U(lenU,M2),VT(lenVT,N2),S(N2))
       !Allocate array needed by LAPACK
       allocate(WORK(LWMAX))

       !Arrays used to compute Hessian Vector product
       ! H^{-1}*A^{T}b = V*SI*UT*A^{T}b
       allocate(SI(lenH,N2),VSI(N2,N2),VSIUT(N2,lenU),xsol(N2,1))
       allocate(V(N2,lenVT),UT(M2,lenU))
       allocate(ATb(lenU,1))

       !Here I read the Matrix H which is my Hessian
       iunit = 81
       open(iunit,file='hessian_column.bin',status='old',&
  &    action='read',form='unformatted',access='direct',&
  &    recl=m2)
       do j=1,m2
         read(iunit,rec=j) H(:,j)
       enddo
       close(iunit)

     !If option == 1 Compute the SVD
     if (option .eq. 1) then

       !You need first to run a test of memory
       !this first call to SGESDD with LWORK = -1 is the test
       LWORK = -1
       CALL SGESDD( 'Singular vectors', M2, N2, H, lenH, S, U, lenU, VT,&
  &                lenVT, WORK, LWORK, IWORK, INFO )

       !Give the correct size to the work array
       !according to the memory test done above
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       LWORK = INT(WORK(1))

       !Print to check the correct size of array work2
       print *, LWMAX, WORK(1), LWORK
       allocate(work2(LWORK))
!
!     Compute SVD.
      CALL SGESDD( 'Singular vectors', M2, N2, H, lenH, S, U, lenU, VT,&
  &                lenVT, WORK2, LWORK, IWORK, INFO )

        !Transpose of V^{T} and U to get V and U^{T}
        do i=1,N2
         V(:,i) = VT(i,:)
         UT(i,:) = U(:,i)
        enddo

        !Save singular values and singular vectors in a file
        open(iunit,file='singular_values.txt',status='unknown')
        do j=1,m2
          write(iunit,*) S(j)
        enddo
        close(iunit)
        open(22,file='V.matrix',status='new')
        open(33,file='UT.matrix',status='new')
         do i=1,N2
          do j=1,N2
           write(22,*) V(i,j)
           write(33,*) UT(i,j)
          enddo
         enddo
        close(22)
        close(33)

       !If option == 2 do not compute SVD but read singular
       !values and singular vectors from files
       elseif ( option .eq. 2 ) then

       !Read singular values and singular vectors
       open(22,file='V.matrix',status='old')
       open(33,file='UT.matrix',status='old')
        do i=1,N2
         do j=1,N2
          read(22,*) V(i,j)
          read(33,*) UT(i,j)
         enddo
        enddo
       close(22)
       close(33)
       open(iunit,file='singular_values.txt',status='old')
        do j=1,m2
         read(iunit,*) S(j)
        enddo
       close(iunit)

       endif

       !Read vector A^{T}b
       ATb(:,:) = 0.
       open(iunit,file='atb.bin',status='old',&
  &    action='read',form='unformatted',access='direct',&
  &    recl=m2)
       read(iunit,rec=1) ATb(:,1)
       close(iunit)
       !Change sign 
       xsol = ATb
       ATb = xsol * -1.

       !Choose how many singular values you want
       write(*,*) 'Singular values to use values /', m2, ':'
       read(*,*) val      

       !Flush all the arrays
       SI(:,:) = 0.
       VSI(:,:) = 0.
       VSIUT(:,:) = 0.
       xsol(:,1) = 0.

       !Inverse of S
       do i=1,val
        SI(i,i) = 1. / S(i)
       enddo 

        m = N2
        k = lenVT
        n = N2
        !First multiplication
        !V*S^-1
        call sgemm('N','N',m,n,k,alpha,V,m,SI,k,beta,VSI,m)
        !Second multiplication
        !V*S^{-1}*U^{T}
        m = N2
        k = N2
        n = N2
        call sgemm('N','N',m,n,k,alpha,VSI,m,UT,k,beta,VSIUT,m)
        !Third multiplication
        !xsol = solution = H^{-1}*A^{T}b = V*S^{-1}*U^{T}*A^{T}b
        m = N2
        k = N2
        n = 1
        call sgemm('N','N',m,n,k,alpha,VSIUT,m,ATb,k,beta,xsol,m)

        !This is the way you call the multiplication of matrices
        !C = alpha*A*B + beta*C
        !A m x k    B k x n    C m x n
        !CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M ) 

       !Write the vector solution xsol = H^{-1}*ATb
       open(iunit,file='xsol.out',status='unknown',&
  &    action='write')
       do i=1,N2
        write(iunit,*) xsol(i,1)
       enddo
       close(iunit)

       print *, 'End of computation'
       deallocate(SI,VSI,VSIUT)
       deallocate(xsol,ATb)
       deallocate(H,U,VT)
       deallocate(S,WORK)
       if (option .eq. 1) then 
          deallocate (WORK2)
       endif
       deallocate(iwork)
       deallocate(V,UT)
 
       endprogram matrix_svd






      subroutine read_time(tractionvec,nsstk,nsdip,msub,ncomp,stcomp,ntime)

       !COMMON VARIABLES
       IMPLICIT NONE
       integer bit, iunit, i, j, k, l, reclent
       integer,intent(inout) :: nsstk, nsdip, msub, stcomp, ncomp, ntime
       real,intent(inout) :: tractionvec(ntime,msub*stcomp*ncomp)
       integer col1

       bit=1
       ! Total elements used for traction matrix (frequency)
       col1=nsstk*nsdip*stcomp*ncomp/4
       reclent=ntime*col1

       iunit=16
        write(*,*) ' Reading Green functions from dat/TRACT_time.bin '
            OPEN(iunit,file='dat/TRACT_time.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
            do i=1,4
            read(iunit,rec=i) tractionvec(:,1+(col1*(i-1)):col1*i)
            enddo
       close(iunit)




!      LWORK = -1
!      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
!  &                LDVT, WORK, LWORK, IWORK, INFO )
!      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!     Compute SVD.
!
!      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT,&
!  &                LDVT, WORK, LWORK, IWORK, INFO )


      end subroutine read_time

