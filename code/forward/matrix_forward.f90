        program forward

    use iso_fortran_env
    implicit none

    integer :: i4
    integer(kind=int8)  :: i8
    integer(kind=int16) :: i16
    integer(kind=int32) :: i32
    integer(kind=int64) :: i64

       integer ntime, ncomp, msub, stcomp
       integer lensyn, lenslip, nx, ny, i, j, k
       real,dimension(:,:),allocatable :: tractionvec, a, source
       real scalfac, dx, dy
       real, dimension(:,:), allocatable :: x, b, AT
       integer iunit
       integer ii, jj, kk, iii, len, m, n
       real al, be

      INTEGER          M2, N2
      INTEGER          LDA2, LDU2, LDVT2
      INTEGER          LWMAX
      integer*8  lt8
      INTEGER          INFO, LWORK
      INTEGER,dimension(:),allocatable :: IWORK
      REAL,DIMENSION(:),ALLOCATABLE :: S2, WORK, WORK2
      REAL,DIMENSION(:,:),ALLOCATABLE :: A2, U2, VT2


      m2 = 34560
      n2 = 34560
      lda2 = m2
      ldu2 = m2
      ldvt2 = n2
      lwmax = 1000
      allocate(IWORK(8*N2))
      allocate(A2(LDA2,N2),U2(LDU2,M2),VT2(LDVT2,N2),S2(N2))
      allocate(WORK(LWMAX))

       write(*,*) 'Time samples in Green functions:'
       !read(*,*) ntime
       write(*,*) 'Subfaults msub:'
       !read(*,*) msub
       write(*,*) 'Length of synthetic seismogram:'
       !read(*,*) lensyn
       write(*,*) 'Length of slip-rate functions:'
       !read(*,*) lenslip
       ntime = 301
       msub = 64
       lensyn = 351
       lenslip = 180

       ncomp = 3
       stcomp = 378*3

       allocate(tractionvec(ntime,msub*ncomp*stcomp))
       allocate(a(lensyn,lenslip*msub*3))
       allocate(b(lensyn,1),x(lenslip*3*msub,1))
       allocate(source(lenslip,msub))
       allocate(at(lenslip*msub*3,lensyn))


  print *, lensyn,lenslip*msub*3
       !Read slip-rate functions as a 2D array
       iunit=20
       open(iunit,file='source.srcf',status='old',action='read')
       do i=1,lenslip
        read(iunit,*) source(i,:)
       enddo
       close(iunit)


       !Arrange the source into a 1D vector
       j=1
       k=1
       do i=1,64
        j=1+(i-1)*180*3
        k=j+179
        x(j:k,1) = source(1:180,i)*-1.
       enddo



       k = 180*64*3
       !do i=1,k
       !  write(34,*) x(i,1)
       !enddo

       !Size of subfaults along strike and dip
       write(*,*) 'dx and dy in grid (km):'
       !read(*,*) dx, dy
       dx = 500. !dx*1000.
       dy = 500. !dy*1000.
       write(*,*) 'nx and ny coarse grid (subfaults):'
       !read(*,*) nx, ny
       nx = 8
       ny = 8
       !Read Green's funcions
       call read_time(tractionvec,nx,ny,msub,ncomp,stcomp,ntime)


       A(:,:) = 0.

       iii=1
       kk=1
       do k=1,64
        do j=1,3
         ii=1
         jj=301
         kk=(k-1)*378*9
         do i=1,180
       !  print *, ii, ii+250, 51, 301
         !p=[ii 301 1 jj]
          A(ii:ii+251,iii)=tractionvec(50:301,kk+j)
          jj=jj-1
          ii=ii+1
          iii=iii+1
         enddo
        enddo
       enddo
       do i=1,301
       write(22,*) A(i,25), tractionvec(i,73)
       enddo

       b(:,1) = 0.
       al = 1.
       be = 0.
       k = 180*64*3
       m = lensyn
       n = 1
       A = A*(dx*dy*0.02) 

       do i=1,lensyn
        AT(:,i) = A(i,:)
       enddo
       !do i=1,lensyn
       ! write(21,*) AT(1,i)
       !enddo


       !Forward modeling     Ax = b
       call sgemm('N','N',m,n,k,al,A,m,x,k,be,b,m)

       A2(:,:) = 0.
       m = 180*64*3
       k = lensyn
       n = 180*64*3
       !Forward modeling     Ax = b
       call sgemm('N','N',m,n,k,al,AT,m,A,k,be,A2,m)

!      CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M ) 
       do i=1,351
         write(34,*) b(i,1)
       enddo
       do i=1,m
         write(35,*) A2(i,1)
       enddo


       LWORK = -1
      CALL SGESDD( 'Singular vectors', M2, N2, A2, LDA2, S2, U2, LDU2, VT2,&
  &                LDVT2, WORK, LWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LWORK = INT(WORK(1))

!      l32 = int(work(1))
print *, LWMAX, WORK(1), LWORK, l32
stop
      allocate(work2(LWORK))
!
!     Compute SVD.
!
!      CALL SGESDD( 'Singular vectors', M2, N2, A2, LDA2, S2, U2, LDU2, VT2,&
!  &                LDVT2, WORK2, LWORK, IWORK, INFO )

       do i=1,m2
         write(37,*) s2(i)
       enddo


        deallocate(A2,U2,VT2)
        deallocate(S2,WORK,WORK2)
        deallocate(iwork)
        deallocate(AT)
        deallocate(b,x)
        deallocate(tractionvec,a,source)
        endprogram forward






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

