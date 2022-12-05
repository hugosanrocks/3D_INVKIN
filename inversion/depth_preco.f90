      !Depth precondition subroutines 
      !for non-progressive inversion

      subroutine depth_preco_go(green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: i, ii, j, k, iunit

      if (green_mesh%rake_opt .eq. 1) then
       !Multiply model by 1/z      
       !z has to be in meters I guess
       k = 1
       do i=1,green_mesh%msub
        do j=1,green_mesh%interp_i
         green_mesh%model1(k) = green_mesh%model1(k) / (green_mesh%fault(i,3)**green_mesh%depth_coef)
         green_mesh%grad1(k) = green_mesh%grad1(k) * (green_mesh%fault(i,3)**green_mesh%depth_coef) 
         k = k + 1
        enddo
       enddo
      elseif (green_mesh%rake_opt .eq. 2) then
       k=1
       do i=1,green_mesh%msub
        do ii=1,2
         do j=1,green_mesh%interp_i
          green_mesh%model2(k) = green_mesh%model2(k) / (green_mesh%fault(i,3)**green_mesh%depth_coef)
          green_mesh%grad2(k) =  green_mesh%grad2(k) * (green_mesh%fault(i,3)**green_mesh%depth_coef)
          k=k+1
         enddo
        enddo
       enddo
      endif

      endsubroutine depth_preco_go


      !Depth precondition subroutines 
      !for non-progressive inversion

      subroutine  depth_preco_back(green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: i, ii, j, k, iunit
      real,dimension(:),allocatable :: vector    ! current gradient


      if (green_mesh%rake_opt .eq. 1) then
       !Multiply model by 1/z      
       !z has to be in meters I guess
       k = 1
       do i=1,green_mesh%msub
        do j=1,green_mesh%interp_i
         green_mesh%model1(k) = green_mesh%model1(k) * (green_mesh%fault(i,3)**green_mesh%depth_coef)
         k = k + 1
        enddo
       enddo
      elseif (green_mesh%rake_opt .eq. 2) then
       k=1
       do i=1,green_mesh%msub
        do ii=1,2
         do j=1,green_mesh%interp_i
          green_mesh%model2(k) = green_mesh%model2(k) * (green_mesh%fault(i,3)**green_mesh%depth_coef)
          k=k+1
         enddo
        enddo
       enddo
      endif

      endsubroutine depth_preco_back



      !Subroutine reading fault geometry

      subroutine read_fault(green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: i, iunit

      !Read geometry of nodes
      iunit=23
      open(iunit,file=green_mesh%out//'fault.mus',status='old',&
  &        action='read')
      do i=1,green_mesh%msub
        read(iunit,*) green_mesh%fault(i,:)
      enddo
      close(iunit)

      endsubroutine read_fault


