      subroutine  interp_axi(proc_mesh)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!=======================================================================================
         IMPLICIT NONE
         INCLUDE 'proc.h'
         TYPE (mesh) :: proc_mesh

!        Variables used in this subroutine
         REAL t_sim, t_slip, y1(6), y2(6), const, slope, corr
         REAL x(2), y(2)
         INTEGER cont, cont2, i, j8, j9

!        Interpolation if simulation_dt is less than slip-rate_dt
         if (proc_mesh%simdt .lt. proc_mesh%slipdt) then

         j8=2        !linear regresion constants j8 j9
         j9=0

         !Set to zero the array not to have trash values
         proc_mesh%stinterp(:,:)=0.d0
 

         !IMPORTANT
         !For time t=0.0 stress values are set to zero and saved to file

         !Cycle over the time simulation
         cont2=1                               ! Counter
         t_sim=proc_mesh%simdt                ! first time sample (simulation)
         DO WHILE (t_sim .le. proc_mesh%simt) ! do untiol final simulation time

         
!         Cycle over the slip rate time
          t_slip=0.d0
          DO cont=1,proc_mesh%slipsam

!         Condition to detect if the stress values for t_slip are
!         near a value of t_sim
          if ((t_slip .le. t_sim+(proc_mesh%simdt/2.d0)).and.&
   &      (t_slip .ge. t_sim-(proc_mesh%simdt/2.d0))) then
 
!         Condition to detect the upper and lower limit of the interpolation
          if (t_slip .le. t_sim) then
!          print *, t_sim-proc_mesh%simdt, t_slip, t_sim, cont2-1, cont, cont2
!          print *, cont2

          y1=proc_mesh%stsi(cont2-1,:)   !read the lower values of stress
          y2=proc_mesh%stsi(cont2,:)     !read the upper values of stress
          x(1)=t_sim-proc_mesh%simdt
          x(2)=t_sim

!         For each value TAU TAU' TAU'' SIGMA.. used the upper and lower
!         values of interpolation as "Y" values and the simul_time as "X"
          do i=1,6
           y(1)=y1(i)
           y(2)=y2(i)

!         perform a linear regresion with two points (x1,y1) (x2,y2)
          call linreg(j8,j9,const,slope,corr,x,y)

!         estimate the interpolated values for TAU TAU' TAU'' SIGMA..          
          proc_mesh%stinterp(cont,i)=((t_slip*slope) + const)
        !  if (cont .eq. 34) then
        !  print *, 'const slope corr x y', const,slope,corr,x,y
        !  print *, cont, x(1), y1(i), x(2), y2(i),t_slip,proc_mesh%stinterp(cont,i)
        !  endif!          print *,t_sim-proc_mesh%simdt, y(1)
!         print *,t_sim, y(2)
!          print *, const, slope, i,  t_slip, interp
          enddo

!         Second condition to detect the up and down limit of the interpolation
          else if (t_slip .ge. t_sim) then
!          print *, t_sim, t_slip, t_sim+proc_mesh%simdt, cont2, cont, cont2+1
!          print *, cont2i
          y1=proc_mesh%stsi(cont2,:)
          y2=proc_mesh%stsi(cont2+1,:)
          x(1)=t_sim
          x(2)=t_sim+proc_mesh%simdt

!         For each value TAU TAU' TAU'' SIGMA.. used the upper and lower
!         values of interpolation as "Y" values and the simul_time as "X"
          do i=1,6
           y(1)=y1(i)
           y(2)=y2(i)

!         perform a linear regresion with two points (x1,y1) (x2,y2)
          call linreg(j8,j9,const,slope,corr,x,y)

!         estimate the interpolated values for TAU TAU' TAU'' SIGMA..          
          proc_mesh%stinterp(cont,i)=((t_slip*slope) + const)
        !  if (cont .eq. 34) then
        !  print *, 'const slope corr x y', const,slope,corr,x,y
        !  print *, cont, x(1), y1(i), x(2), y2(i),t_slip,proc_mesh%stinterp(cont,i)
        !  endif!          print *,t_sim-proc_mesh%simdt, y(1)
!          print *,t_sim, y(2)
!          print *, const, slope, i,  t_slip, interp
          enddo


          endif
          
          !number of interpolated values, number of time samples in the range
          proc_mesh%interp_i=cont


          endif

          t_slip = t_slip + proc_mesh%slipdt  ! Check for the next slip time step
          
          ENDDO

         t_sim = t_sim + proc_mesh%simdt      ! Check for the next simulation time step
         cont2 = cont2 + 1
         ENDDO


         else
         print *, 'Error slip-rate dt > simulation dt'


         endif





         end subroutine interp_axi
