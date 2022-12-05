         subroutine arrange_grad(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem(1,2), i, k, n1, n2

         !Flush the array
         green_mesh%slipr2(:,:) = 0.d0


         mem(1,1) = 1
         mem(1,2) = 1+green_mesh%interp_i

         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%Interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            green_mesh%slipr2(k,i) = green_mesh%grad2(n1)
            green_mesh%slipr2(k,i+green_mesh%msub) = green_mesh%grad2(n2)
          enddo
         enddo

         if (green_mesh%lap_opt .eq. 2) then
          call smooth_grad(green_mesh)
         elseif (green_mesh%lap_opt .eq. 1) then
          call smooth_lap(green_mesh)
         endif

         !Return additional gradient to 1D array (interp_i*msub*2)
         !=======================================!
         green_mesh%gradad(:) = 0.d0
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem(1,1) + green_mesh%Interp_i*2*(i-1)+(k-1)
            n2 = mem(1,2) + green_mesh%interp_i*2*(i-1)+(k-1)
            green_mesh%grad2(n1) = green_mesh%slipr2(k,i)
            green_mesh%grad2(n2) = green_mesh%slipr2(k,i+green_mesh%msub)
          enddo
         enddo

         endsubroutine arrange_grad







         subroutine arrange_grad1d(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem, i, k, n1

         !Flush the array
         green_mesh%slipr2(:,:) = 0.d0
         mem = 1

         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
            n1 = mem + green_mesh%Interp_i*(i-1)+(k-1)
            green_mesh%slipr2(k,i) = green_mesh%grad1(n1)
          enddo
         enddo


         if (green_mesh%lap_opt .eq. 2) then
          call smooth_grad1d(green_mesh)
         elseif (green_mesh%lap_opt .eq. 1) then
          call smooth_lap(green_mesh)
         endif

         !Return additional gradient to 1D array (interp_i*msub*2)
         !=======================================!
         green_mesh%gradad(:) = 0.d0
         do k=1,green_mesh%interp_i
          do i=1,green_mesh%msub
           n1 = mem + green_mesh%Interp_i*(i-1)+(k-1)
           green_mesh%grad1(n1) = green_mesh%slipr2(k,i)
          enddo
         enddo


         endsubroutine arrange_grad1d




         subroutine smooth_grad(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer m, n, i, j, k, ii, add_i
         real, dimension(:,:), allocatable :: matrixin, matrixout
         integer r1, r2, i1, i2, n1, n2, repeat

         !print *, 'smoothing gradient'

         allocate(matrixin(green_mesh%nsdip,green_mesh%nsstk))
         allocate(matrixout(green_mesh%nsdip,green_mesh%nsstk))
         matrixin(:,:) = 0.
         matrixout(:,:) = 0.
         r1 = green_mesh%r1        !radius on dip direction
         r2 = green_mesh%r2        !radius on strike direction
         repeat = green_mesh%repeat    !times repeating the smoothing filter
         n1 = green_mesh%nsdip
         n2 = green_mesh%nsstk

         do ii=1,2      ! 2D components
          add_i = (ii-1)*green_mesh%msub
          !parts of slip vector
          do j=1,green_mesh%interp_i
           k=1
           do m = 1,green_mesh%nsdip
            do n = 1,green_mesh%nsstk
             matrixin(m,n) = green_mesh%slipr2(j,k+add_i)
             k = k+1
            enddo
           enddo

            !Convolution with triangle functions
            do i=1,repeat
             do i2=1,n2
              call triangle(r1,n1,matrixin(:,i2),matrixout(:,i2))
             enddo
             do i1=1,n1
              call triangle(r2,n2,matrixout(i1,:),matrixin(i1,:))
             enddo
            enddo
            k=1
            do m = 1,green_mesh%nsdip
             do n = 1,green_mesh%nsstk
              green_mesh%slipr2(j,k+add_i) = matrixin(m,n)
              k = k+1
             enddo
            enddo

          enddo
         enddo


         deallocate(matrixin,matrixout)
         endsubroutine smooth_grad






         subroutine smooth_grad1d(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer m, n, i, j, k
         real, dimension(:,:), allocatable :: matrixin, matrixout
         integer r1, r2, i1, i2, n1, n2, repeat

         !print *, 'smoothing gradient 1D'

         allocate(matrixin(green_mesh%nsdip,green_mesh%nsstk))
         allocate(matrixout(green_mesh%nsdip,green_mesh%nsstk))
         matrixin(:,:) = 0.
         matrixout(:,:) = 0.
         r1 = green_mesh%r1        !radius on dip direction
         r2 = green_mesh%r2        !radius on strike direction
         repeat = green_mesh%repeat    !times repeating the smoothing filter
         n1 = green_mesh%nsdip
         n2 = green_mesh%nsstk

          !parts of slip vector
          !print *, m, n, i
          do j=1,green_mesh%interp_i
           k=1
           do m = 1,green_mesh%nsdip
            do n = 1,green_mesh%nsstk
             matrixin(m,n) = green_mesh%slipr2(j,k)
             !write(54,*) green_mesh%slipr2(j,k)
             k = k+1
            enddo
           enddo

            !Convolution with triangle functions
            do i=1,repeat
             do i2=1,n2
              call triangle(r1,n1,matrixin(:,i2),matrixout(:,i2))
             enddo
             do i1=1,n1
              call triangle(r2,n2,matrixout(i1,:),matrixin(i1,:))
             enddo
            enddo
            k=1
            do m = 1,green_mesh%nsdip
             do n = 1,green_mesh%nsstk
              green_mesh%slipr2(j,k) = matrixin(m,n)
              write(62,*) green_mesh%slipr2(j,k)
              k = k+1
             enddo
            enddo
          enddo

         deallocate(matrixin,matrixout)
         endsubroutine smooth_grad1d




!==================================================================
subroutine boxconv(nbox, nx, xx, yy)
  implicit none
  integer,intent(in):: nx,nbox
  integer::i,ny
  real,dimension(nx),intent(in)::xx
  real,dimension(nx+nbox-1),intent(out)::yy
  real,dimension(:),allocatable::bb

  allocate(bb(nx+nbox))
  if (nbox < 1 .or. nbox > nx) then
     write(0,*) "boxconv: error in the length of input!"
  endif
  ny = nx+nbox-1
  do i= 1, ny
     bb(i) = 0.
  enddo
  bb(1) = xx(1)
  do i= 2, nx
     bb(i) = bb(i-1) + xx(i)  ! make B(Z) = X(Z)/(1-Z)
  enddo
  do i= nx+1, ny
     bb(i) = bb(i-1)
  enddo
  do i= 1, nbox
     yy(i) = bb(i)
  enddo
  do i= nbox+1, ny
     yy(i) = bb(i) - bb(i-nbox) ! make Y(Z) = B(Z)*(1-Z**nbox)
  enddo
  do i= 1, ny
     yy(i) = yy(i) / nbox
  enddo
  deallocate(bb)
end subroutine boxconv

!=================================================================
!apply triangle filter to input xx: yy=triangle*xx
subroutine triangle(nbox,nd,xx,yy)
  implicit none

  integer, intent(in)::nbox,nd
  integer::i,np,nq
  real,dimension(nd),intent(in)::xx
  real,dimension(nd),intent(out)::yy
  real,dimension(:),allocatable::pp,qq

  allocate(pp(nd+nbox-1))
  allocate(qq(nd+2*nbox-2))
  call boxconv(nbox,nd,xx,pp)
  np=nbox+nd-1
  call boxconv(nbox,np,pp,qq)
  nq=nbox+np-1
  do i=1,nd
     yy(i)=qq(i+nbox-1)
  enddo
  do i=1,nbox-1
     yy(i)=yy(i)+qq(nbox-i) !fold back near end
  enddo
  do i=1,nbox-1
     yy(nd-i+1)=yy(nd-i+1)+qq(nd+(nbox-1)+i) !fold back far end
  enddo
  deallocate(pp)
  deallocate(qq)
end subroutine triangle


         subroutine smooth_lap(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh


         !call lap_filter(green_mesh)
         call lap_filter_2d(green_mesh)

         endsubroutine smooth_lap




