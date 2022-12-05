       subroutine model_d(x,y,lx,ly,matrix,dir_n,detect)

       !Projection of a 3D vector (x,y,z)
       !onto a plane according to strike and dip
       !vectors (matrix)

       implicit none
       !Dimension of vectors
       integer, intent(inout) :: lx, ly, dir_n
       !Vectors to be convolved
       real, intent(inout) :: x(lx*ly*3), y(lx*ly*2), matrix(dir_n,2,3)
       integer, intent(inout) :: detect(ly)
       real vec(3), vec2(2)
       !Counters
       integer i, k, cont, cont2


       do i=1,ly
        do k=1,lx
         vec(1)=x( k    + ((lx*3)*(i-1))  )
         vec(2)=x( k+lx + ((lx*3)*(i-1))  )
         vec(3)=x(k+lx*2+ ((lx*3)*(i-1))  )
         !print *, k    + ((lx*3)*(i-1)), k+lx + ((lx*3)*(i-1)), k+lx*2+ ((lx*3)*(i-1))
         !Projection along strike and dip
         vec2(:)=0.
         do cont=1,2
          do cont2=1,3
           vec2(cont)= vec2(cont)+ &
  &                   matrix(detect(i),cont,cont2)*vec(cont2)
          enddo
         enddo
         y(k    + ((lx*2)*(i-1)) ) = vec2( 1  )
         y(k+lx + ((lx*2)*(i-1)) ) = vec2( 2  )
         !print *, k+ ((lx*2)*(i-1)), k+lx +((lx*2)*(i-1))
        enddo
       enddo
    


       end subroutine model_d



       subroutine model_c(x,y,lx,ly,matrix,dir_n,detect)

       !Projection of a 2D vector, along strike and dip
       !to a 3D (x,y,z) vector

       implicit none
       !Dimension of vectors
       integer, intent(inout) :: lx, ly, dir_n
       !Vectors to be convolved
       real, intent(inout) :: x(lx*ly*2), y(lx*ly*3), matrix(dir_n,2,3)
       integer, intent(inout) :: detect(ly)
       real vec(3), vec2(2)
       !Counters
       integer i, k, cont, cont2

       
       do i=1,ly
        do k=1,lx
         vec2(1)=x( k    + ((lx*2)*(i-1))  )
         vec2(2)=x( k+lx + ((lx*2)*(i-1))  )
         !print *, k+ ((lx*2)*(i-1)), k+lx +((lx*2)*(i-1))
         !Projection along strike and dip
         vec(:)=0.
         do cont=1,3
          do cont2=1,2
           vec(cont)=vec(cont)+ &
  &                  matrix(detect(i),cont2,cont)*vec2(cont2)
          enddo
         enddo
         y(k     + ((lx*3)*(i-1)) ) = vec( 1 )
         y(k+lx  + ((lx*3)*(i-1)) ) = vec( 2 )
         y(k+lx*2+ ((lx*3)*(i-1)) ) = vec( 3 )
         !print *, k+ ((lx*3)*(i-1)), k+lx +((lx*3)*(i-1)), k+lx*2+ ((lx*3)*(i-1))
        enddo
       enddo
        

       end subroutine model_c

