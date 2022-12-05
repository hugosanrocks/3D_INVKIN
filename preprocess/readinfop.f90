      subroutine read_infop(proc_mesh)

         IMPLICIT NONE
         INCLUDE 'proc.h'
         TYPE (mesh) :: proc_mesh

         integer iunit, i
        
         !Directory containing input data files P_C* TUA_C* SIGMA_**_C*
         proc_mesh%dat='dat/'
         !Directory containing output data files
         proc_mesh%out='out/'

         iunit=10

         !Read reciever's geometry
         print *, 'Reading simul.info'
         open(iunit,file=proc_mesh%dat//'simul.info',status='old',action='read')
         read(iunit,*) proc_mesh%rt, proc_mesh%ot
         read(iunit,*) proc_mesh%nsta, proc_mesh%ncomp, proc_mesh%msub
         read(iunit,*) proc_mesh%simsam, proc_mesh%simt, proc_mesh%simdt
         read(iunit,*) proc_mesh%stress_opt
         print *, 'Simulation dt', proc_mesh%simdt
         close(iunit)

         !Read focal mechanism information 
         open(iunit,file=proc_mesh%dat//'focal.info',status='old',action='read')
         read(iunit,*) proc_mesh%dir_n
         allocate(proc_mesh%vnorm(proc_mesh%dir_n,3))
         allocate(proc_mesh%vslip(proc_mesh%dir_n,3))
         do i = 1, proc_mesh%dir_n
         proc_mesh%dir_i = i
         read(iunit,*) proc_mesh%stk, proc_mesh%dip, proc_mesh%rak
         print *, 'Focal mechanism', i, proc_mesh%stk, proc_mesh%dip, proc_mesh%rak
         call norm_vec(proc_mesh)
         enddo
         close(iunit)


         !Read slip rate point source information
         open(iunit,file=proc_mesh%dat//'sliprate.info',status='old',action='read')
         read(iunit,*) proc_mesh%slipsam, proc_mesh%slipt
         close(iunit)
         proc_mesh%slipdt=proc_mesh%slipt/(proc_mesh%slipsam - 1)
         print *, 'Slip rate samples and slipdt', proc_mesh%slipsam, proc_mesh%slipdt

         open(iunit,file=proc_mesh%dat//'filpreproc.info',&
  &           status='old',action='read')
         read(iunit,*) proc_mesh%f_pre
         read(iunit,*) proc_mesh%fc_pre
         read(iunit,*) proc_mesh%order
         close(iunit)


      end subroutine read_infop
