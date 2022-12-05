       subroutine edge(green_mesh)

        !COMMON VARIABLES
        IMPLICIT NONE
        INCLUDE 'green.h'
        TYPE (mesh) :: green_mesh

        integer i, iunit
        real, dimension(:), allocatable :: vector_in

!        real dist(green_mesh%msub,green_mesh%msub)
!        real lambda, midc(1,3), x(1,24), y(1,12), z(1,12)
!        real expx(1,24), expy(1,12), expp(12,24), al, be, mat(288,288)
!        integer nx1, ny1, nx, ny, nxy, mm


!        mm=green_mesh%msub
!        al=1.
!        be=0.
        !SAME AS Mathilde Radiguet THESIS EQ (3.3) 
        !Controling parameters of covariance matrix
!        lambda = 1500.     !distance of exponential decay

!        green_mesh%fault(:,:)=0.
!        dist(:,:) = 0.

         !File where to read the subfault positions
!         iunit=10
!         open(iunit,file=green_mesh%dat//'fault.pos',status='!old',&
! &            action='read',access='DIRECT',recl=green_mesh%msub*4*green_mesh%ncomp)
!         !Read subfault positions (x,y,z)
!         read(iunit,rec=1) green_mesh%fault(:,:)
         !print *, green_mesh%fault(1,:)
!         close(iunit)
!         nx1 = 4
!         ny1 = 4
!         nx = nx1 / 2
!         ny = ny1 / 2
!         nxy = nx1*(ny-1) + nx
!         midc(1,1) = (green_mesh%fault(nxy,1) + green_mesh%fault(nxy+1,1)) / 2.
!         midc(1,2) = (green_mesh%fault(nxy,2) + green_mesh%fault(nxy+nx1,2)) / 2.
!         midc(1,3) = (green_mesh%fault(nxy,3) + green_mesh%fault(nxy+nx1,3)) / 2.
         !print *, midc, 'middle point'

         !checar valores aqui
!         do i=1,nx1
!          x(1,i) = midc(1,1) -midc(1,1) - real(nx)*1.5 +0.75 + real(i-1)*1.5
!         enddo
!         do i=1,ny1
!          y(1,i) = midc(1,2) -midc(1,2) - real(ny)*1.5 + 0.75 + real(i-1)*1.5
!         enddo
!         x(1,:) = x(1,:) * !1000.
!         y(1,:) = y(1,:) * 1000.

!         do i=1,nx1
!          expx(1,i)= 1. / exp(-1.*x(1,i) / lambda )
!          expx(1,i)= expx(1,i) + 1. / exp( x(1,i)/ lambda );
!         enddo
!         expx(:,:) = expx(:,:) / maxval(expx(1,:))

!         do i=1,ny1
!          expy(1,i)= 1. / exp(-1.*y(1,i) / lambda )
!          expy(1,i)= expy(1,i) + 1. / exp( y(1,i)/ lambda );
!         enddo
!         expy(:,:) = expy(:,:) / maxval(expy(1,:))

       !Estimate distance and correlation matrices
!        do j=1,nx1
!         expp(:,j) = expy(1,:)
!        enddo
!       do i=1,ny1
!         expp(i,:) = expp(i,:) + expx(1,:)
!       enddo
 
!       k=1
!       do i=1,ny1
!        do j=1,nx1
!         green_mesh%ce(k,1) = expp(i,j)
!         k = k + 1
!        enddo
!       enddo
!       do i=1,green_mesh%msub
!         green_mesh%ce(:,i) = green_mesh%ce(:,1)
!       enddo

       allocate(vector_in(green_mesh%msub))

       vector_in(:) = 0.
       green_mesh%ce(:,:) = 0.
       iunit = 111
       open(iunit,file='dat/mat.dat',status='unknown')
       do i=1,green_mesh%msub
        read(iunit,*) vector_in(1:green_mesh%msub)
        green_mesh%ce(i,:) = vector_in(:)
       enddo

!            call sgemm('N','T',mm,mm,mm,al,        &
!      &       green_mesh%ce,mm,green_mesh%ce,mm,be,mat,mm)

       !green_mesh%ce(:,:) = mat(:,:)/mat(1,1)
       !do j=1,green_mesh%msub
       !do i=1,green_mesh%msub
       !write(33,*) green_mesh%ce(i,j)
       !enddo
       !enddo
       
       endsubroutine edge
