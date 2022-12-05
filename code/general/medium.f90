        subroutine id_mu(proc_mesh)

        !COMMON VARIABLES
         IMPLICIT NONE
         INCLUDE 'proc.h'
         TYPE (mesh) :: proc_mesh

         !Variables needed by this subroutine
         INTEGER i, j, nlayer, model_opt, iunit
         REAL,DIMENSION(:,:),ALLOCATABLE :: medium
         REAL,DIMENSION(:),ALLOCATABLE :: mu, lambda
         integer,dimension(:,:),allocatable :: me


         allocate(me(proc_mesh%msub,6))

         !Flush arrays
         proc_mesh%mus(:)= 0.
         proc_mesh%fault(:,:) = 0.

         !File where to read the subfault positions
         iunit=10
         open(iunit,file=proc_mesh%dat//'fault.pos',status='old',&
 &            action='read',access='DIRECT',recl=proc_mesh%msub*4*proc_mesh%ncomp)
         !Read subfault positions (x,y,z)
         read(iunit,rec=1) proc_mesh%fault(:,:)
         !print *, proc_mesh%fault(:,:)
         close(iunit)


         !Read medium properties (Several options)
         open(iunit,file=proc_mesh%dat//'medium.info',status='old',&
 &            action='read')
         read(iunit,*) model_opt, nlayer
         allocate(medium(nlayer,6),mu(nlayer),lambda(nlayer))
         !Model option 1 = Homegeneous medium
         if (model_opt .eq. 1) then 
          !Read subfault positions (z,Vp,Vs,rho,Qp,Qs)
          read(iunit,*) medium(:,:)
            !Estimate Lame parameters at each layer	
            mu(nlayer) = ((medium(nlayer,3)*1000.)**2)*(medium(nlayer,4)*1000.)
            lambda(nlayer) = ((medium(nlayer,2)*1000.)**2)*(medium(nlayer,4)*1000.) - (2. * (mu(nlayer)*1000.))
         elseif (model_opt .eq. 2) then
          do i=1,nlayer
           read(iunit,*) medium(i,:)
            mu(i) = ((medium(i,3)*1000.)**2)*(medium(i,4)*1000.)
            lambda(i) = ((medium(i,2)*1000.)**2)*(medium(i,4)*1000.) - (2. * mu(i))
            write(*,*) ' Layer:', i, 'mu:', mu(i)
          enddo
         elseif (model_opt .eq. 3) then
           print *, 'Heterogeneous media, not implemented yet'
         else
           print *, 'Wrong option to read velocity model'
         endif
         close(iunit)

        ! print *, mu(:), 'muss'
        ! print *, model_opt, proc_mesh%fault
         if (model_opt .eq. 1) then
         proc_mesh%mus(:) = mu(nlayer)
         elseif (model_opt .eq. 2) then
          do i=1,nlayer-1
           do j=1,proc_mesh%msub
             if ((proc_mesh%fault(j,3) .ge. medium(i,1)) &
 &              .and. (proc_mesh%fault(j,3) .lt. medium(i+1,1))) then
                proc_mesh%mus(j) = mu(i)

                me(j,1) = i

             elseif (proc_mesh%fault(j,3) .ge. medium(nlayer,1) ) then
                proc_mesh%mus(j) = mu(nlayer)

                me(j,1) = nlayer

             elseif ( proc_mesh%fault(j,3) .lt. 0. ) then
                print *, 'Wrong subfault location (negative depth). Check fault.pos'
             endif
           enddo
          enddo
         endif


         open(iunit,file=proc_mesh%out//'fault.mus',status='unknown')
         do i=1,proc_mesh%msub
           write(iunit,*) proc_mesh%fault(i,:), me(i,1) !proc_mesh%mus(i)
         enddo
         close(iunit)


         deallocate(me)

         deallocate(medium,mu,lambda)
         endsubroutine id_mu
