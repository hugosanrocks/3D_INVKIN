SUBROUTINE progressbar(i)
integer, intent(inout) :: i
! Uses FLUSH (f2003 but not f95) and assumes unit 6 = stdout
CALL sleep(1)
WRITE(6,"(A)",ADVANCE='NO')'*'
FLUSH 6
END SUBROUTINE progressbar
