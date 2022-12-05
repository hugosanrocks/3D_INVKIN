          program newread

       implicit none
       integer i,j,k,n,m,nsta,nsub,ncyc,nt,iunit,iunit2
       REAL, DIMENSION(:,:), ALLOCATABLE :: dat
       REAL, DIMENSION(:), ALLOCATABLE :: datvec
       INTEGER, DIMENSION(:,:), ALLOCATABLE :: index

       character*4, DIMENSION(:), ALLOCATABLE :: file1
       character*2 num


          iunit=10
          iunit2=11


          nsub=2       !subfaults
          nsta=5       !receivers
          nt=9670      !time steps

          allocate(dat(nt,nsub),index(nsta,nsub),datvec(nt*nsub))

          !Pointer to know where to write the data read
          DO i = 1,nsta
            k = 0
            DO j = 1,nsub
                Index(i,j) = i + (k*nsta)
                k = k + 1
            ENDDO
          ENDDO

!======================== first component

!======================== P

          open(iunit2,file='P_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
!             open(iunit,file='Px'//num//'',status='old',&
!   &         form='unformatted',access='direct',recl=nt*nsub*4)

              read(iunit,rec=1) datvec

!             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               !write(iunit2,rec=index(n,j)) dat(:,j)
               print *, index(n,j)
               enddo




          enddo

          close(iunit2)




!======================== component

!======================== TAU_P

          open(iunit2,file='TAU_P_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
             open(iunit,file='TAU_Px'//num//'',status='old',&
   &         form='unformatted',access='direct',recl=nt*nsub*4)

             read(iunit,rec=1) datvec

             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               write(iunit2,rec=index(n,j)) dat(:,j)
               enddo




          enddo

          close(iunit2)







!======================== component

!======================== TAU_PP

          open(iunit2,file='TAU_PP_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
             open(iunit,file='TAU_PPx'//num//'',status='old',&
   &         form='unformatted',access='direct',recl=nt*nsub*4)

             read(iunit,rec=1) datvec

             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               write(iunit2,rec=index(n,j)) dat(:,j)
               enddo




          enddo

          close(iunit2)







!======================== component

!======================== SIGMA_XY

          open(iunit2,file='SIGMA_XY_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
             open(iunit,file='SIGMA_XYx'//num//'',status='old',&
   &         form='unformatted',access='direct',recl=nt*nsub*4)

             read(iunit,rec=1) datvec

             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               write(iunit2,rec=index(n,j)) dat(:,j)
               enddo




          enddo

          close(iunit2)







!======================== component

!======================== SIGMA_XZ

          open(iunit2,file='SIGMA_XZ_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
             open(iunit,file='SIGMA_XZx'//num//'',status='old',&
   &         form='unformatted',access='direct',recl=nt*nsub*4)

             read(iunit,rec=1) datvec

             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               write(iunit2,rec=index(n,j)) dat(:,j)
               enddo




          enddo

          close(iunit2)






!======================== component

!======================== SIGMA_YZ

          open(iunit2,file='SIGMA_YZ_C1',status='unknown',&
   &      form='unformatted',access='direct',recl=nt*4)


          do n = 1,nsta

             write(num,'(I2.2)') n
             open(iunit,file='SIGMA_YZx'//num//'',status='old',&
   &         form='unformatted',access='direct',recl=nt*nsub*4)

             read(iunit,rec=1) datvec

             close(iunit)

             k = 1
             do i = 1,nt
               do m = 1,nsub
                dat(i,m) = datvec(k)
                k = k + 1
               enddo
             enddo


               do j = 1,nsub
               write(iunit2,rec=index(n,j)) dat(:,j)
               enddo




          enddo

          close(iunit2)





          deallocate(dat,datvec,index)
          endprogram newread
