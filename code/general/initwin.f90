      subroutine initwin(green_mesh)

!----------------------Definition of variables--------------------------------------
      !COMMON VARIABLES
      implicit none
      !Define all the variables needed to read stresses
      ! and calculate the tractions associated "Green's
      ! functions", variables in include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      real :: time4samples


      !if iregular mesh is not used
      if ( green_mesh%opt_mesh .lt. 3 ) then

      !if it is used
      elseif ( green_mesh%opt_mesh .eq. 3 ) then
         green_mesh%msub = green_mesh%nnodes
      endif


      !1D vectors used for optimization toolbox
      !for progressive inversion
      if (green_mesh%for_opt .eq. 2) then
       !time samples in source functions
       green_mesh%interp_i = green_mesh%samp4win(green_mesh%synwin) + green_mesh%mext
       write(6, *) ' Time samples (interp_i):', green_mesh%interp_i
       time4samples = real(green_mesh%interp_i-1)*green_mesh%slipdt
       green_mesh%ntint = int(time4samples/green_mesh%intdt)+1
      !if inversion is not PIS
      else
       !Nothing interp_i and ntint come from simul.info
       write(6, *) ' Time samples (interp_i):', green_mesh%interp_i
      endif

      !green_mesh%ntint = green_mesh%ntint + green_mesh%mxint
      write(6, *) ' Time samples for interp (ntint):', green_mesh%ntint

      !Model sizes
      ! according to a 3D vector (x,y,z)
      green_mesh%modelsize = green_mesh%interp_i * green_mesh%ncomp * &
  &                          green_mesh%msub
      green_mesh%sizeint =  green_mesh%ntint * green_mesh%ncomp * &
  &                         green_mesh%msub

      !according to a 2D vector (strike,dip)
      green_mesh%modelsize2 = green_mesh%interp_i*2*green_mesh%msub
      !according to a 1D vector (rake fixed)
      green_mesh%modelsize1 = green_mesh%interp_i*green_mesh%msub

      !Size for convolution array
      green_mesh%lensyn=green_mesh%ntint+green_mesh%trac_i-1

      !Array to store traction on the fault (gradient)
      allocate(green_mesh%tottrac(green_mesh%ntint,green_mesh%ncomp*green_mesh%msub))
      call monitor_alloc('27%tottrac',' allo')

      !Array to store synthetic seismograms
      allocate(green_mesh%syn(green_mesh%lenobs,green_mesh%stcomp))
      call monitor_alloc('28%syn    ',' allo')

      !2 slipmod = matrix containing slip module for each subfualt
      allocate(green_mesh%slipmod(green_mesh%interp_i,green_mesh%msub))
      call monitor_alloc('29%slipmod',' allo')

      !3 4 5
      allocate(green_mesh%model(green_mesh%modelsize))
      call monitor_alloc('30%model3d',' allo')
      allocate(green_mesh%modelint(green_mesh%sizeint))
      call monitor_alloc('31%modelin',' allo')
      allocate(green_mesh%model2(green_mesh%modelsize2))
      call monitor_alloc('32%model2d',' allo')
      allocate(green_mesh%grad2(green_mesh%modelsize2))
      call monitor_alloc('33%grad2d ',' allo')
      allocate(green_mesh%model1(green_mesh%modelsize1))
      call monitor_alloc('34%model1d',' allo')
      allocate(green_mesh%grad1(green_mesh%modelsize1))
      call monitor_alloc('35%grad1d ',' allo')
      allocate(green_mesh%gradad1(green_mesh%modelsize1))
      call monitor_alloc('36%gradad1',' allo')

      !if rake is fixed
      if (green_mesh%rake_opt .eq. 1) then
        allocate(green_mesh%target(green_mesh%interp_i*green_mesh%msub))
        call monitor_alloc('37%target1',' allo')
      elseif (green_mesh%rake_opt .eq. 2) then
        allocate(green_mesh%target(green_mesh%interp_i*green_mesh%msub*2))
        call monitor_alloc('38%target2',' allo')
      endif

      allocate(green_mesh%tracint(green_mesh%interp_i,green_mesh%ncomp*green_mesh%msub))
      call monitor_alloc('39%tracint',' allo')

      !arrays used for matrix multiplication in the regularization
      allocate(green_mesh%slipr(green_mesh%interp_i,green_mesh%ncomp))
      call monitor_alloc('40%slipr  ',' allo')

      !7 Matrix with slip-rate time history in 3D for all subfaults
      allocate(green_mesh%slip(green_mesh%ntint,green_mesh%ncomp*green_mesh%msub))
      call monitor_alloc('41%slip   ',' allo')


      allocate(green_mesh%gradad(green_mesh%modelsize2))
      call monitor_alloc('42%gradad2',' allo')

      allocate(green_mesh%slipr2(green_mesh%interp_i,green_mesh%msub*2))
      call monitor_alloc('43%slipr2 ',' allo')

      allocate(green_mesh%modelp(green_mesh%modelsize))
      call monitor_alloc('44%modelp ',' allo')

      !22 A priori model in 2D
      !allocate(green_mesh%model2p(green_mesh%modelsize2))

      allocate(green_mesh%diag(green_mesh%msub*green_mesh%interp_i))
      call monitor_alloc('45%diag   ',' allo')
      allocate(green_mesh%rup(green_mesh%msub*green_mesh%interp_i))
      call monitor_alloc('46%rup    ',' allo')
      allocate(green_mesh%q_up(green_mesh%modelsize1))
      call monitor_alloc('47%q_up   ',' allo')

      allocate(green_mesh%filter(green_mesh%nsstk,green_mesh%nsdip,green_mesh%interp_i))
      call monitor_alloc('48%filter ',' allo')
      allocate(green_mesh%resfil(green_mesh%nsstk,green_mesh%nsdip,green_mesh%interp_i))
      call monitor_alloc('49%resfil ',' allo')
      allocate(green_mesh%rake_p(green_mesh%modelsize2))
      call monitor_alloc('50%rake_p ',' allo')

      allocate(green_mesh%p_model1d(green_mesh%interp_i,green_mesh%msub))
      call monitor_alloc('51%pmodel1',' allo')
      allocate(green_mesh%p_model2d(green_mesh%interp_i,green_mesh%msub*2))
      call monitor_alloc('52%pmodel2',' allo')


      allocate(green_mesh%values(green_mesh%lensyn,green_mesh%msub))
      call monitor_alloc('53%values ',' allo')



      end subroutine initwin




      subroutine destroywin(green_mesh)

!----------------------Definition of variables--------------------------------------
      !COMMON VARIABLES
      implicit none
      !Define all the variables needed to read stresses
      ! and calculate the tractions associated "Green's
      ! functions", variables in include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh


      deallocate(green_mesh%tottrac)
      call monitor_alloc('27%tottrac',' deal')
      deallocate(green_mesh%syn)
      call monitor_alloc('28%syn    ',' deal')
      deallocate(green_mesh%slipmod)
      call monitor_alloc('29%slipmod',' deal')
      deallocate(green_mesh%model)
      call monitor_alloc('30%model3d',' deal')
      deallocate(green_mesh%modelint)
      call monitor_alloc('31%modelin',' deal')
      deallocate(green_mesh%model2)
      call monitor_alloc('32%model2d',' deal')
      deallocate(green_mesh%grad2)
      call monitor_alloc('33%grad2d ',' deal')
      deallocate(green_mesh%model1)
      call monitor_alloc('34%model1d',' deal')
      deallocate(green_mesh%grad1)
      call monitor_alloc('35%grad1  ',' deal')
      deallocate(green_mesh%gradad1)
      call monitor_alloc('36%gradad1',' deal')
      deallocate(green_mesh%target)
      call monitor_alloc('37/8target',' deal')
      deallocate(green_mesh%tracint)
      call monitor_alloc('39%tracint',' deal')
      deallocate(green_mesh%slipr)
      call monitor_alloc('40%slipr  ',' deal')
      deallocate(green_mesh%slip)
      call monitor_alloc('41%slip   ',' deal')
      deallocate(green_mesh%gradad)
      call monitor_alloc('42%gradad2',' deal')
      deallocate(green_mesh%slipr2)
      call monitor_alloc('43%slipr2 ',' deal')
      deallocate(green_mesh%modelp)
      call monitor_alloc('44%modelp ',' deal')
      deallocate(green_mesh%diag)
      call monitor_alloc('45%diag   ',' deal')
      deallocate(green_mesh%rup)
      call monitor_alloc('46%rup    ',' deal')
      deallocate(green_mesh%q_up)
      call monitor_alloc('47%q_up   ',' deal')
      deallocate(green_mesh%filter)
      call monitor_alloc('48%filter ',' deal')
      deallocate(green_mesh%resfil)
      call monitor_alloc('49%resfil ',' deal')
      deallocate(green_mesh%rake_p)
      call monitor_alloc('50%rake_p ',' deal')
      deallocate(green_mesh%p_model1d)
      call monitor_alloc('51%pmodel1',' deal')
      deallocate(green_mesh%p_model2d)
      call monitor_alloc('52%pmodel2',' deal')
      deallocate(green_mesh%values)
      call monitor_alloc('53%values ',' deal')

      endsubroutine destroywin
