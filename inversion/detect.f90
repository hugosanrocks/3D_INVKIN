     subroutine detect(green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: n, iunit, i, j, k, l, cont                                        ! dimension of the problem
      integer :: predsize
      real,dimension(:),allocatable :: grad_preco, model, grad     ! current gradient
      character*4 :: FLAG                                          ! communication FLAG 


       if (green_mesh%rake_opt .eq. 1) then
         call read_modelpri1d(green_mesh)
       elseif (green_mesh%rake_opt .eq. 2) then
         call read_modelpri(green_mesh)
       endif
       print *, ' Started model'
       !Number of windows growing the model space
       !call prediction(green_mesh)

       !call localmax(green_mesh)


       endsubroutine detect



     subroutine localmax(green_mesh)

      implicit none
!     Define all the variables needed to read models and
!     associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

!     Variables needed only here
      integer :: n, iunit, i, j, k, l, cont                                        ! dimension of the problem
      integer :: predsize
      real,dimension(:),allocatable :: model, maxva
      integer,dimension(:),allocatable :: idmax

        !Read form file and identify the subfaults where
        !the prediction will be done
        cont = 0
        j = 0
        do i=1,green_mesh%msub
         if (green_mesh%flag_d .eq. 'PRED') then
!        if (green_mesh%win(i,green_mesh%dowin) .eq. green_mesh%map(green_mesh%dowin) ) then
          if ( (green_mesh%win(i,1) .ge. green_mesh%dowin-4) &
  &          .and. ( green_mesh%win(i,1) .le. green_mesh%dowin) .and. &
                     (green_mesh%synwin .lt. 8)  )then
             cont = cont + 1
             j = j + 1
             green_mesh%idsub(j) = i
!print *, 'cond 1', i
          elseif ( (green_mesh%win(i,2) .gt. 1 ) & !green_mesh%dowin-5) &
  &          .and. ( green_mesh%win(i,1) .le. green_mesh%dowin-1) .and. &
                     (green_mesh%synwin .eq. 8)  )then
             cont = cont + 1
             j = j + 1
             green_mesh%idsub(j) = i
!print *, 'cond 2', i
          else
          endif
         elseif (green_mesh%flag_d .eq. 'DETC') then
          cont = cont+1
         endif
        enddo
print *, cont, 'cont'

                 !subf * comp * samples
        predsize =  cont
print *, predsize, 'predsize'

        !Array to store the norm of the slip vector
        allocate(model(green_mesh%interp_i),maxva(predsize),idmax(predsize))
        model(:) = 0.
        maxva(:) = 0.
        idmax(:) = 0
 
      !Arrange slip rate model in 1D vector
      !Remove subfaults not to update
      do i=1,green_mesh%msub
          !Select only subfaults near hypocenter
        if (green_mesh%flag_d .eq. 'PRED') then
  
        elseif (green_mesh%flag_d .eq. 'DETC') then
          cont=1+(i-1)*2*green_mesh%interp_i
          !Remove subfaults not to be updated
          model(:) = 0.
          do j=1,2
            do k=1,green_mesh%interp_i
              model(k) = model(k) + green_mesh%model2(cont)**2
              cont = cont + 1
            enddo
          enddo
          !Take the norm 
          model(:) = sqrt(model)
          maxva(i) = maxval(model)
          do n = 1,green_mesh%interp_i
           if ( maxva(i) .eq. model(n)) then
             idmax(i) = n
           endif
          enddo
        endif
      print *, i, idmax(i), maxva(i)
      enddo

      n = maxval(idmax)
      print *, 'n', n
      do i=1,predsize
           !Substract 30 samples
           if (idmax(i) .lt. 30) then
             idmax(i) = 1
           else
             if (n .le. 20) then
              idmax(i) = idmax(i) - n
             elseif (n .gt. 20) then
               idmax(i) = idmax(i) - 20
             endif
           endif
           green_mesh%rsamp(i) = idmax(i)
      print *, i, green_mesh%rsamp(i), 'changed'
      enddo

      do i=1,green_mesh%msub
        write(32,*) green_mesh%rsamp(i)
      enddo      
      call time_mask(green_mesh)

      deallocate(model)
      deallocate(maxva)
      deallocate(idmax)
     endsubroutine localmax

