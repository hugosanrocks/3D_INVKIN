       subroutine stress_axi(proc_mesh)

       !READING STRESS INPUT DATA FROM AXITTRA

       IMPLICIT NONE
       INCLUDE 'proc.h'
       TYPE (mesh) :: proc_mesh

       integer*8 n, reclen
       integer iunit2

        proc_mesh%stsi(:,:) = 0.d0   !Flush the stress input array

        !proc_mesh%file_s='stress_C'//proc_mesh%comp//'.txt'

        proc_mesh%file_s='stress.bin'

        !Open the file to be read stress_C*  *= Component of non-zero force
        iunit2=15

         ! number of record inside stress bin file
         n = proc_mesh%comp_i + (proc_mesh%sta_i-1)*proc_mesh%ncomp

          open(iunit2,file=proc_mesh%dat//proc_mesh%file_s,&
   &      status='old',form='unformatted',access='direct',recl=1024*4*6)

          !Read the stress data from axitra
          read(iunit2,rec=n) proc_mesh%stsi
          print *, n, 'record'
          !Print only to check the read is correct
          !do i=1,proc_mesh%simsam
          !  write(22,*) proc_mesh%stsi(i,:)
          !enddo

       endsubroutine stress_axi
