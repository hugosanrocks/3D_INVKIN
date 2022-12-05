      subroutine npow2(n,npo2)

      integer, intent(inout) :: n, npo2

      npo2 = 2
      DO WHILE (npo2 .lt. n) ! do until next pow2 is found
       npo2 = npo2*2
      ENDDO

      end subroutine npow2
