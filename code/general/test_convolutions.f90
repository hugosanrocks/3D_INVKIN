      program test1

      use omp_lib
      implicit none
      integer i, j, k, n, td_id, chunk, map(2,5)
      real array1(5,2), array2(20,10), array3(24,2)

      chunk = 1

      array1(:,1) = [0., 1., 2., 1., 0.]
      array1(:,2) = [0., 0., 1., 0., 0.]
      do i = 1, 20
       do j = 1, 10
        array2(i,j) = real(i)*real(j)
       enddo
      enddo

      k = 1
      do i = 1, 2
       do j = 1, 5
        map(i,j) = k
        k = k + 1
       enddo
      enddo

      array3(:,:) = 0.
      !$OMP PARALLEL SHARED(array2,array1,array3,map,chunk) PRIVATE(td_id,i)
       td_id = OMP_get_thread_num()
       print *, 'Thread starting: ', td_id
       call convolutions(map,array1,array2,array3,td_id)
      !$OMP END PARALLEL

print *, array3(:,1)

      endprogram test1

      subroutine convolutions(map,array1,array2,array3,td_id)

      implicit none
      integer i, j, k, nsta, nsub, lenx, leny, lenz
      real, intent(inout) :: array1(5,2), array2(20,10), array3(24,2)
      integer, intent(inout) :: map(2,5), td_id
      real x(5), y(20), z(24), zz(24)

      zz(:) = 0.
      nsta = 2
      nsub = 5
      lenx = 5
      leny = 20
      lenz = 24

      zz(:) = 0.
      !$OMP DO REDUCTION(+:array3) SCHEDULE(STATIC,1)
       do i = 1, nsta
         x(:) = array1(:,i)
        do j = 1, nsub
         y(:) = array2(:,map(i,j))
         !IF ((i .eq. 2) .and. (j .lt. 3)) then
          !print *, 'y', map(i,j), y
          !print *, i, map(i,j), td_id
         !endif
         call conv(x,lenx,y,leny,z,lenz)
         zz(:) = zz(:) + z(:)
         !IF ((i .eq. 1) .and. (j .eq. 2)) then
         ! print *, 'z',z
         !endif
        enddo
        array3(:,i) = zz
       enddo
      !$OMP END DO

      end subroutine
