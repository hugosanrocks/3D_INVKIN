

         subroutine mul_src(proc_mesh)

  USE INIT_BUTTERWORTH_MOD
  USE FILTFILT_BUTTERWORTH_MOD

         IMPLICIT NONE
         INCLUDE 'proc.h'
         TYPE (mesh) :: proc_mesh

         integer*8 n, reclen, nfile
         integer iunit2, icont, mfault
         real, dimension(:), allocatable :: vecread
         !------------------------------------------
         !     variable declaration for butterworth
         !-------------------------------------------
           TYPE(butter) :: butt
           INTEGER nstep
           REAL dt
           INTEGER nrec
         !-----------------------------

         !Option to read stress data from GEODG3D

         proc_mesh%stsi(:,:) = 0.d0   !Flush the stress input array
         n=proc_mesh%simsam           !Number of samples in Green's functions
         allocate(vecread(n))         !Memony needed for reading records
         !Jump of samples inside the files
         icont=proc_mesh%mjump/(n)
         mfault=(icont/(proc_mesh%nsta))+1       !Number of fault used to identify
                                                 !medium properties (mu)
         !Byte length of each stress input element TAU TAU' TAU'' SIGMA..
         INQUIRE(iolength=reclen) vecread(1)

         do nfile=1,6  ! six files to read per source-station
           !choose the file name P_SXXX_CX or TAU_P_SXXX_CX, etc.
           if (nfile .eq. 1) then
             proc_mesh%file_s='P_C'//proc_mesh%comp//''
           elseif(nfile .eq. 2) then
             proc_mesh%file_s='TAU_P_C'//proc_mesh%comp//''
           elseif(nfile .eq. 3) then
             proc_mesh%file_s='TAU_PP_C'//proc_mesh%comp//''
           elseif(nfile .eq. 4) then
             proc_mesh%file_s='SIGMA_XY_C'//proc_mesh%comp//''
           elseif(nfile .eq. 5) then
             proc_mesh%file_s='SIGMA_XZ_C'//proc_mesh%comp//''
           elseif(nfile .eq. 6) then
             proc_mesh%file_s='SIGMA_YZ_C'//proc_mesh%comp//''
           endif
           !Open the file to be read TAU TAU' TAU'' SIGMA..
           iunit2=15
           OPEN(iunit2,FILE=proc_mesh%dat//proc_mesh%file_s,&
    &      FORM='UNFORMATTED',STATUS='UNKNOWN',ACCESS='DIRECT',&
    &      recl=reclen*n)
           !Read the stress component for each subfault, station, component
           read(iunit2,rec=proc_mesh%sta_i+icont) vecread
           close(iunit2)
           !-------------------------------------------------------------------
           !OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
           !-------------------------------------------------------------------
           !Parameters needed
           butt%order = 2
           butt%fc = 0.5 !4.
           nstep = proc_mesh%simsam
           nrec  = proc_mesh%stcomp
           dt    = proc_mesh%simdt
           !Asign the time series to be filtered
           proc_mesh%tseries(:) = vecread(:)
           !Initialize the butterworth filter
           CALL INIT_BUTTERWORTH(butt, dt)
           ! Apply the butterworth filter
           CALL FILTFILT_BUTTERWORTH(proc_mesh%tseries, butt, nstep)
           vecread(:) = proc_mesh%tseries(:)
           ! Print the resultsprint input and filtered signals to output file
           !!OPEN(UNIT=1,FILE='test_butterworth.txt',STATUS='UNKNOWN')
           !!t=tbegin
           !!DO ii = 1,nstep
           !! WRITE(1,'(E12.5,x,E12.5,x,E12.5)') t, tseries_total(ii), copy_tseries(ii,2)
           !! t=t+dt
           !!ENDDO
           !!CLOSE(unit=1)
           !!PRINT *,' '
           !!PRINT *,'Results in file test_butterworth.txt :)'
           !-----------------------------------------------------------------
           ! END OF OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
           !-----------------------------------------------------------------

          !Save filtered strees component
          proc_mesh%stsi(:,nfile) = vecread(:)
          !print *, 'rec=', proc_mesh%sta_i+icont
         enddo        !END CYCLE OVER 6 FILES TO BE READ

       !Scale the stress tensor
       !Divide the stress tensor by the corresponding mu
       proc_mesh%stsi(:,:)=proc_mesh%stsi(:,:)/proc_mesh%mus(mfault)


       !Deallocate the required memory
       deallocate(vecread)
       endsubroutine mul_src



