      subroutine detect_focmec(proc_mesh)

      implicit none
      !Define all the variables needed to read stresses and calculate the tractions
      !associated "Green's functions", variables in include/green.h
      INCLUDE 'proc.h'
      TYPE (mesh) :: proc_mesh

      !Variables needed only by this program
      integer i, iunit

      iunit=25
      open(iunit,file=proc_mesh%dat//'/map_focal.dat',status='old',&
  &        action='read')

      do i=1,proc_mesh%msub
        read(iunit,*) proc_mesh%vnorm_i(i) 
      enddo
      close(iunit)

      endsubroutine detect_focmec
