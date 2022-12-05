      subroutine ortogonal(green_mesh)

!----------------------Definition of variables--------------------------------------

      implicit none

!     Define all the variables needed to read stresses and calculate the tractions
!     associated "Green's functions", variables in include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: iunit2, i, j, k              ! counters
      real,dimension(:,:),allocatable :: x    ! current state matrices
      real,dimension(:),allocatable :: y      ! current state matrices
      REAL pi/3.1415926535897932/

      iunit2=34


!         green_mesh%stk=green_mesh%stk*(2.0*pi/360.0)
!         green_mesh%dip=green_mesh%dip*(2.0*pi/360.0)
!         green_mesh%rak=green_mesh%rak*(2.0*pi/360.0)


         !Unitary slip vector used for transformation
         !(Stein p.218 eq.1)
!         green_mesh%vslip(1) =-1.0*(cos(green_mesh%rak)*sin(green_mesh%stk)) + &   !-1.0
!  &      sin(green_mesh%rak)*cos(green_mesh%dip)*cos(green_mesh%stk)
!         green_mesh%vslip(2) =cos(green_mesh%rak)*cos(green_mesh%stk) + &    !stk 0 -----> 180
!  &      sin(green_mesh%rak)*cos(green_mesh%dip)*sin(green_mesh%stk)
!         !Because of the orientation on the mesh
!         green_mesh%vslip(3) = sin(green_mesh%rak)*sin(green_mesh%dip)  !-1 z para arriba
      !Read focal mechanism information 
!!      open(iunit,file=green_mesh%dat//'focal.info',&
!!  &        status='old',action='read')
!!       read(iunit,*) green_mesh%stk, green_mesh%dip, green_mesh%rak
!!      close(iunit)

!!         call norm_vec(green_mesh)
         !print *, 'normal vector',green_mesh%vnorm


        !Save memory to read the model and gradient info
        allocate(x(green_mesh%interp_i*green_mesh%msub,green_mesh%ncomp),&
    &            y(green_mesh%interp_i*green_mesh%msub))

        OPEN(iunit2,FILE=green_mesh%dat//'slip_xyz.ascii',&
    &   status='unknown')


        !Read model and gradient info as matrices
          do i=1,green_mesh%msub
           do j=1,green_mesh%interp_i
            read(iunit2,*) x(i,:)
            y(i)= x(i,1)*green_mesh%vnorm(green_mesh%vnorm_i(i),1) + &
  &               x(i,2)*green_mesh%vnorm(green_mesh%vnorm_i(i),2) + &
  &               x(i,3)*green_mesh%vnorm(green_mesh%vnorm_i(i),3)
           enddo
          enddo
          close(iunit2)

        OPEN(iunit2,FILE=green_mesh%dat//'ortogonal.ascii',&
    &   status='unknown')

        !Read model and gradient info as matrices
          do i=1,green_mesh%interp_i*green_mesh%msub
           write(iunit2,*) y(i)
          enddo
        close(iunit2)



         deallocate(x,y)
         end subroutine ortogonal
