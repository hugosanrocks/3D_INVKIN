
         subroutine vectors_time_filter1_in(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer i, n1, n2
         real tseries_in(green_mesh%interp_i)
         
          n1 = (green_mesh%cont_i-1)*green_mesh%interp_i !first sample
          n2 = n1 + green_mesh%interp_i                  !last sample
          green_mesh%tseries_in(:) = green_mesh%grad1(n1:n2)

          endsubroutine vectors_time_filter1_in


         subroutine vectors_time_filter1_out(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer i, n1, n2


          !Return additional gradient to 1D array (interp_i*msub*2)
          !=======================================!
            n1 = (green_mesh%cont_i-1)*green_mesh%Interp_i
            n2 = n1 + green_mesh%interp_i
            green_mesh%grad1(n1:n2) = green_mesh%tseries_out(:)


         endsubroutine vectors_time_filter1_out


      subroutine time_filter_preco(green_mesh)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      integer :: i, j, n1, n2, operation_type, ii, jj
      real h, tol
      real,dimension(:,:),allocatable :: vector_in_2d, vector_out_2d
      REAL,DIMENSION(:,:),ALLOCATABLE :: lx,lz,sigma,theta

      n1 = green_mesh%nsdip
      n2 = green_mesh%nsstk
      green_mesh%h_lap = green_mesh%stk_s
      h = green_mesh%h_lap
      operation_type = 1

!      call read_vectors_lap2d(n1,n2,green_mesh%lx,green_mesh%lz,green_mesh%sigma,&
!  &                           green_mesh%theta,green_mesh%tol_lap)

      !To test use an spike
      !green_mesh%vector_in(:,:) = 0.
      !green_mesh%vector_in(4,7) = 0.
      !green_mesh%vector_in(9,4) = 1.

      do i=1,green_mesh%interp_i

      !Write to check the gradient before filtering
      !do ii=1,n1
      !  do jj=1,n2
      !   write(777,*) green_mesh%vector_in(ii,jj)
      !  enddo
      !enddo

       green_mesh%cont_i = i
       call vectors_2d_smooth1_in(green_mesh)

       call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &        green_mesh%sigma,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)
       !open(11,file='output_vector_2d_oncethrough',access='direct',status='replace',recl=4*n1*n2)
       !write(11,rec=1)vector_out_2d(:,:)
       !close(11)
       green_mesh%vector_in(:,:)= green_mesh%vector_out(:,:)
       call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &         green_mesh%sigma,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)

       call vectors_2d_smooth1_out(green_mesh)

      !Write only to check the filter effect
!      do ii=1,n1
!        do jj=1,n2
!         write(666,*) green_mesh%vector_out(ii,jj)
!        enddo
!      enddo

      enddo


      endsubroutine time_filter_preco


      subroutine init_filter_time_preco(tc1,tc2,dt,n,filter)

      implicit none
      integer, intent(inout) :: n
      real, intent(inout) :: tc1, tc2, dt
      real, intent(out) :: filter(n)
      real t(n)
      integer i, j

      dt = 0.3
      !sam=50; n
      tc1=0.5;
      tc2=3;
      t(:) = 0.
      do i=2,n-1
      t(i) = t(i-1)+dt
       if (t(i-1) <= tc1) then
          filter(i-1) = exp(-((t(i-1)-tc1)**2/tc1**2))
       elseif (t(i-1) > tc1) then
          filter(i-1) = exp(-((t(i-1)-tc1)**2/tc2**2))
       endif
      enddo

      do i=1,n
       write(667,*) filter(i)
      enddo

      endsubroutine init_filter_time_preco
