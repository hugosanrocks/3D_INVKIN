            program readas

            real stsi
            integer iunit, iunit2, re, i, n, nsta, m, k
            character (len=1) s
            nsta=5
            n=6447
            m=n*nsta

            re=4
            iunit=10
            iunit2=11
       OPEN(iunit2,FILE='P',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='P'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)

            print *, '***********************************************'
            do j=1,n

             read(iunit,rec=j) stsi
             print *, k, stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)


       OPEN(iunit2,FILE='TAU_P',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='TAU_P'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            do j=1,n

             read(iunit,rec=j) stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)



       OPEN(iunit2,FILE='TAU_PP',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='TAU_PP'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            do j=1,n

             read(iunit,rec=j) stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)




       OPEN(iunit2,FILE='SIGMA_XY',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='SIGMA_XY'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            do j=1,n

             read(iunit,rec=j) stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)





       OPEN(iunit2,FILE='SIGMA_XZ',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='SIGMA_XZ'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            do j=1,n

             read(iunit,rec=j) stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)





       OPEN(iunit2,FILE='SIGMA_YZ',&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            k=1
            do i=1,nsta

            write(s,'(I1.1)') i
       OPEN(iunit,FILE='SIGMA_YZ'//s,&
    &  FORM='UNFORMATTED',&
    &  STATUS='UNKNOWN',ACCESS='DIRECT',recl=re)


            do j=1,n

             read(iunit,rec=j) stsi

             write(iunit2,rec=k) stsi
            k=k+1
            enddo

           close(iunit)

            enddo

           close(iunit2)








           end program readas
