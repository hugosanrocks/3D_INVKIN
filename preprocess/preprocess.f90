      subroutine preprocess(proc_mesh)

      implicit none

      !Define all the variables needed to read stresses and calculate the tractions
      !associated "Green's functions", variables in include/green.h
      INCLUDE 'proc.h'
      TYPE (mesh) :: proc_mesh

      !Variables needed only by this program
      integer i, k, m
      real :: start, stop


      WRITE(6 ,*) ''
      WRITE(6, *) '================================================='
      WRITE(6, *) ' PREPROCESS, PREPARING UNITARY TRACTION VECTORS  '
      WRITE(6, *) '================================================='

      !Read reciever's geometry
      call read_infop(proc_mesh)


      !Reading stress tensors
      WRITE(*,*) 'Reading stress input files '

      !6 files with simsam = samples in simulation
      allocate(proc_mesh%stsi(proc_mesh%simsam,6))
      allocate(proc_mesh%stinterp(proc_mesh%simsam,6))
      allocate(proc_mesh%tractionv(proc_mesh%simsam,3))

      !Subfaults's positions array
      allocate(proc_mesh%fault(proc_mesh%msub,3))
      allocate(proc_mesh%mus(proc_mesh%msub))

      !Time series array used to filtered Green's functions
      allocate(proc_mesh%tseries(proc_mesh%simsam))

      !Assign medium properties to each subfault (mu, lambda)
      call id_mu(proc_mesh)

      proc_mesh%tfft_i = 1

      !Time mark
      call cpu_time(start)

      !Read stress files from different codes
      if (proc_mesh%stress_opt .eq. 1) then
      write(6,*) 'Input stress data comes from GEODG3D'
      elseif (proc_mesh%stress_opt .eq. 2) then
      write(6,*) 'Input stress data comes from AXITRA'
      elseif (proc_mesh%stress_opt .eq. 3) then
      write(6,*) 'Input stress data comes from DWN-IBEM'
      else
      write(6,*) 'Wrong inputfile option 1=GEODG3D 2=AXITRA 3=DWN-IBEM'
      endif

      !Allocate array saving index to detect focal mechanism
      allocate(proc_mesh%vnorm_i(proc_mesh%msub))

      write(*,*) ' Detecting focal mechanism from dat/map_focal.info'
      call detect_focmec(proc_mesh)

      !estimate unitary traction vectors
      do m=1,proc_mesh%msub          !number of subfaults
       print *, 'Subfault no: ', m 
       proc_mesh%subf_i = m
       proc_mesh%mjump=(m-1)*proc_mesh%simsam*proc_mesh%nsta   
       !jump inside input files
         proc_mesh%sta_i=1
         do i=1,proc_mesh%nsta
          proc_mesh%comp_i=1
          do k=1,proc_mesh%ncomp
           write(proc_mesh%sta,'(I3.3)') i
           write(proc_mesh%comp,'(I1.1)') k
             if (proc_mesh%stress_opt .eq. 1) then
              !read files
              call mul_src(proc_mesh)
              !Interpolate input files every slip dt
              call interpolation(proc_mesh)     
             elseif (proc_mesh%stress_opt .eq. 2) then !From AXITRA
              !read files
              call stress_axi(proc_mesh)
              !interpolate
              call interp_axi(proc_mesh)
             elseif (proc_mesh%stress_opt .eq. 3) then !From DWN-IBEM
              !read files
              call stress_dwn(proc_mesh)
              !interpolate
              call interp_dwn(proc_mesh)
             else
              print *, 'Wrong option to read stress input files'
             endif
           !4) Construction of stress tensor and
           !estimation of tractions, they are save into TRACT_time.bin
           call traction(proc_mesh)
           proc_mesh%tfft_i = proc_mesh%tfft_i + 1
           proc_mesh%comp_i=proc_mesh%comp_i+1
          enddo
          proc_mesh%sta_i=proc_mesh%sta_i+1
         enddo
      enddo
      call cpu_time(stop)
      write(6,*) '*****************************************'
      write(6,*) ' Time used for preprocess: ', stop-start, 'sec'
      write(6,*) '*****************************************'
      WRITE(6, *) '================================================='
      WRITE(6, *) ' PREPROCESS FINISHED, OUTPUT "dat/TRACT_time.bin"'
      WRITE(6, *) '================================================='

     deallocate(proc_mesh%vnorm_i)
     deallocate(proc_mesh%vnorm,proc_mesh%vslip)
     deallocate(proc_mesh%fault,proc_mesh%mus)
     deallocate(proc_mesh%tractionv)
     deallocate(proc_mesh%stsi,proc_mesh%stinterp)
     deallocate(proc_mesh%tseries)
     end subroutine preprocess
