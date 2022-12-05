          program random_numbers

      !Reset the random seed
      USE init_random_seed_mod

          implicit none
          integer :: lx
          real, dimension(:), allocatable :: x
          integer :: i
          real :: val, a, b

          !how many numbers you want
          lx = 100
          allocate(x(lx))

          !normal dist. parameters
          a = 0.
          b = 1.

          !Where to save values outside the library
          val = 0.

          !reset random seed
          call INIT_RANDOM_SEED()

          !get the random numbers according to normal
          !distribution N(a,b) -> a=mean, b=sigma
          do i=1,lx
           call NORMAL_DIST (a, b, val)
           x(i) = val
          enddo

          deallocate(x)
          endprogram random_numbers
