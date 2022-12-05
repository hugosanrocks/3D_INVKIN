          subroutine random_vector(a,b,x,lx)

      !Reset the random seed
      USE init_random_seed_mod

          implicit none
          integer, intent(inout) :: lx
          real, intent(inout) :: a, b, x(lx,1)
          integer :: i
          real :: val

          !Where to save values
          val = 0.

          call INIT_RANDOM_SEED()

          do i=1,lx
           call NORMAL_DIST (a, b, val)
           x(i,1) = val
          enddo


          endsubroutine random_vector
