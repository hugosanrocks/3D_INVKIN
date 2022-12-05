
      subroutine read_obs2(green_mesh)

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      
       real obs(875,120), o(875)
       integer i, j, k, omp_get_thread_num, ind(3,40), iunit
       character*3 sta
       character*1 comp
       character*21 fil
       character*12 fil2

       do i=1,40
        ind(1,i) = (i-1)*3 +1
        ind(2,i) = (i-1)*3 +2
        ind(3,i) = (i-1)*3 +3
       enddo

!$OMP PARALLEL PRIVATE(sta,comp,iunit,fil,o)
!$OMP DO
       do i=1,40
         do j=1,3
          write(sta,'(I3.3)') i
          write(comp,'(I1.1)') j 
          iunit = omp_get_thread_num()
           !print *, i, sta, iunit, omp_get_thread_num()
           fil='dat/obs_S'//sta//'_C'//comp//'.ascii'
           open(iunit,file=fil,status='old')
             read(iunit,*) o(:)
             green_mesh%obs(:,ind(j,i)) = o(:)
           close(iunit)
         enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL
       print *, green_mesh%obs(170,1:5)

!NOT NEEDED

!!$OMP PARALLEL PRIVATE(sta,comp,iunit)
!!$OMP DO
!!       do i=1,40
!!         write(sta,'(I3.3)') i
!!           iunit = omp_get_thread_num()
!!           print *, i, sta, iunit, fil2, omp_get_thread_num()
!!           fil2='dat/obs_S'//sta
!!           open(iunit,file=fil2,status='unknown')
!!           do k=1,875
!!             write(iunit,*) obs(k,ind(1,i)), obs(k,ind(2,i)), obs(k,(ind(3,i)))
!!           enddo
!!           close(iunit)
!!       enddo
!!$OMP END DO
!!$OMP END PARALLEL

       call weight_cov(green_mesh)

      end subroutine read_obs2
