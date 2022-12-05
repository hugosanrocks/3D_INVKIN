       subroutine monitor_alloc(var,flag)

        integer :: iunit
        character*5 :: flag
        character*10 :: var

        iunit = 22
        open(iunit,file='alloc_monitor.out',&
  &          status='old',access='append')
         write(iunit,*) var, flag
        close(iunit)

       endsubroutine monitor_alloc
