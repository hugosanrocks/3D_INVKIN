

         subroutine grad_2d_smooth(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         !Variables needed only here
         integer :: i, j, i1, i2, n1, n2, ii, repeat
         real h, tol, r1, r2
         real,dimension(:,:),allocatable :: vector_in_2d, vector_out_2d
         REAL,DIMENSION(:,:),ALLOCATABLE :: lx,lz,sigma,theta

         n1 = green_mesh%nsdip
         n2 = green_mesh%nsstk
         r1 = 1.
         r2 = 1.
         repeat = 2

         green_mesh%grad1(:) = 0.
         green_mesh%grad1(29*4+10) = 1.
         do i=1,green_mesh%interp_i*290
          write(54,*) green_mesh%grad1(i)
         enddo


         if (green_mesh%rake_opt .eq. 1) then
          do ii=1,green_mesh%interp_i
           green_mesh%cont_i = ii
           call vectors_2d_smooth1_in(green_mesh)
           do i=1,repeat
            do i2=1,n2
             call triangle(r1,n1,green_mesh%vector_in(:,i2),&
  &               green_mesh%vector_out(:,i2))
            enddo
            do i1=1,n1
             call triangle(r2,n2,green_mesh%vector_out(i1,:),&
  &               green_mesh%vector_in(i1,:))
            enddo
           enddo
           call vectors_2d_smooth1_out(green_mesh)
          enddo
         endif

         do i=1,green_mesh%interp_i*290
          write(55,*) green_mesh%grad1(i)
         enddo

         endsubroutine grad_2d_smooth


         subroutine vectors_2d_smooth1_in(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem, i, j, k, l, cont, n1
         real vec2
         real tc1, tc2, filter(50), dt
         integer n

         mem = 1

          cont = 1
          do i=1,green_mesh%nsdip
           do j=1,green_mesh%nsstk
            n1 = mem + green_mesh%Interp_i*(cont-1)+(green_mesh%cont_i-1)
            green_mesh%vector_in(i,j) = green_mesh%grad1(n1)
            cont = cont + 1
           enddo
          enddo

          endsubroutine vectors_2d_smooth1_in


         subroutine vectors_2d_smooth1_out(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem, i, j, k, l, cont, n1
         real vec2

         mem = 1

          !Return additional gradient to 1D array (interp_i*msub*2)
          !=======================================!
          cont = 1
          do i=1,green_mesh%nsdip
           do j=1,green_mesh%nsstk
            n1 = mem + green_mesh%Interp_i*(cont-1)+(green_mesh%cont_i-1)
            green_mesh%grad1(n1) = green_mesh%vector_out(i,j)
            cont  = cont + 1
           enddo
          enddo


         endsubroutine vectors_2d_smooth1_out




         subroutine vectors_2d_smooth2_in(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem, i, j, k, l, cont, n1
         real vec2
         real tc1, tc2, filter(50), dt
         integer n

         if ( green_mesh%comp_i .eq. 1 ) then
          mem = 1
         elseif ( green_mesh%comp_i .eq. 2 ) then
          mem = 1+green_mesh%interp_i
         endif
! 
        !do k=1,green_mesh%interp_i
          cont = 1
          do i=1,green_mesh%nsdip
           do j=1,green_mesh%nsstk
            n1 = mem + green_mesh%Interp_i*2*(cont-1)+(green_mesh%cont_i-1)
            green_mesh%vector_in(i,j) = green_mesh%grad2(n1)
            cont = cont + 1
           enddo
          enddo

!      tc1 = 0.5
!      tc2 = 3.
!      dt = 0.3
!      n = 50
!      filter(:) = 0.

!      call init_filter_time_preco(tc1,tc2,dt,n,filter)


          endsubroutine vectors_2d_smooth2_in


         subroutine vectors_2d_smooth2_out(green_mesh)

         IMPLICIT NONE
         INCLUDE 'green.h'
         TYPE(mesh) :: green_mesh

         integer mem, i, j, k, l, cont, n1
         real vec2

         if ( green_mesh%comp_i .eq. 1 ) then
          mem = 1
         elseif ( green_mesh%comp_i .eq. 2 ) then
          mem = 1+green_mesh%interp_i
         endif

          !Return additional gradient to 1D array (interp_i*msub*2)
          !=======================================!
          cont = 1
          do i=1,green_mesh%nsdip
           do j=1,green_mesh%nsstk
            n1 = mem + green_mesh%Interp_i*2*(cont-1)+(green_mesh%cont_i-1)
            green_mesh%grad2(n1) = green_mesh%vector_out(i,j)
            cont  = cont + 1
           enddo
          enddo
!         enddo

         endsubroutine vectors_2d_smooth2_out





      subroutine smooth_preco(green_mesh)

      implicit none
      !Define all the variables needed to read models and
      !associated gradient, variables in ../include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !Variables needed only here
      integer :: i, j, n1, n2, operation_type, ii, jj
      real h, tol
      real,dimension(:,:),allocatable :: vector_in_2d, vector_out_2d
      REAL,DIMENSION(:,:),ALLOCATABLE :: lx,lz,sigma,theta

      n1 = green_mesh%nsdip
      n2 = green_mesh%nsstk
      green_mesh%h_lap = green_mesh%stk_s
      h = green_mesh%h_lap
      operation_type = 1

      if (green_mesh%rake_opt .eq. 1) then
      do i=1,green_mesh%interp_i
       green_mesh%cont_i = i
       call vectors_2d_smooth1_in(green_mesh)
       call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &        green_mesh%phi,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)
       green_mesh%vector_in(:,:)= green_mesh%vector_out(:,:)
       call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &         green_mesh%phi,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)
       call vectors_2d_smooth1_out(green_mesh)
      enddo
      elseif (green_mesh%rake_opt .eq. 2) then
      do j=1,2
       green_mesh%comp_i = j
       do i=1,green_mesh%interp_i
        green_mesh%cont_i = i
        call vectors_2d_smooth2_in(green_mesh)
        call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &         green_mesh%phi,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)
        green_mesh%vector_in(:,:)= green_mesh%vector_out(:,:)
        call lap_2d_sparse(green_mesh%vector_in,n1,n2,h,operation_type,green_mesh%lx,green_mesh%lz,&
  &         green_mesh%phi,green_mesh%theta,green_mesh%tol_lap,green_mesh%vector_out)
        call vectors_2d_smooth2_out(green_mesh)
       enddo
      enddo
      else
        write(*,*) 'Error rake option'
        stop
      endif



      endsubroutine smooth_preco


      subroutine read_vectors_lap2d(n1,n2,lx,lz,phi,theta,tol)

IMPLICIT NONE

INTEGER,INTENT(INOUT) :: n1, n2
REAL,INTENT(INOUT) :: lx(n1,n2),lz(n1,n2),phi(n1,n2),theta(n1,n2),tol
INTEGER:: lx_filesize, lz_filesize, theta_filesize, sigma_filesize, use_variable_sigma
CHARACTER(80)::lx_file,lz_file,sigma_file,theta_file

!Read Laplacian Covariance input file
OPEN(72,file='covariance_input')
READ(72,*)lx_file,lz_file,theta_file !HORIZONTAL AND VERTICAL CORRELATION LENGTH VECTOR FILENAMES AND THETA FILE
!!WRITE(*,*)trim(lx_file),trim(lz_file),trim(theta_file)
READ(72,*)tol ! Sparse linear solver tolerance
!!WRITE(*,*)tol
READ(72,*)use_variable_sigma !Use variable sigma or unitary operator (L only)
if (use_variable_sigma==1) then
READ(72,*)sigma_file
end if
CLOSE(72)


write(*,*)" **CHECKING FILE SIZES SMOOTH PRECO** "
!CHECK FILE SIZE OF INPUT FILES TO ENSURE CONSISTENT
INQUIRE(FILE=lx_file, SIZE=lx_filesize)
INQUIRE(FILE=lz_file, SIZE=lz_filesize)
INQUIRE(FILE=theta_file, SIZE=theta_filesize)
lx_filesize=lx_filesize/4
lz_filesize=lz_filesize/4
theta_filesize=theta_filesize/4

if (lx_filesize/=n1*n2) then
print *, 'n1', n1, 'n2', n2
write(*,*)"LX FILE SIZE ERROR"
STOP
end if
if (lz_filesize/=n1*n2) then
write(*,*)"LZ FILE SIZE ERROR"
STOP
end if

 write(*,*)" READING BINARY FILES SMOOTH PRECO "

open(72,file=lx_file,access='direct',recl=4*n1*n2)
read(72,rec=1)lx(:,:)
close(72)

open(72,file=lz_file,access='direct',recl=4*n1*n2)
read(72,rec=1)lz(:,:)
close(72)

open(72,file=theta_file,access='direct',recl=4*n1*n2)
read(72,rec=1)theta(:,:)
close(72)

if (use_variable_sigma==1) then   
   INQUIRE(FILE=sigma_file, SIZE=sigma_filesize)
   if (sigma_filesize/=n1*n2*4) then
      write(*,*)"SIGMA FILE SIZE ERROR"
      STOP
   end if
   open(72,file=sigma_file,access='direct',recl=4*n1*n2)
!   read(72,rec=1)sigma(:,:)
!sigma(:,:) = vector_in(:,:)
   close(72)
end if

open(72,file='sigma_tesst',access='direct',status='replace',recl=4*n1*n2)
write(72,rec=1)phi(:,:)
close(72)

print *, ' *** End reading lap files*** '

      endsubroutine read_vectors_lap2d

subroutine  lap_2d_sparse(vector_in,n1,n2,h,operation_type,lx,lz,sigma,theta,tol,vector_out)

!Modified from Paul Wellington lap_2d_sparse_rot
!OPERATION TYPE
!1-> C_LAP_2D * v (SPARSE SOLVER)
!3-> CINV_LAP_2D * v

IMPLICIT NONE

INTEGER,INTENT(INOUT) :: n1, n2
REAL,INTENT(INOUT) :: vector_in(n1,n2)
REAL,INTENT(INOUT) :: vector_out(n1,n2)
REAL,INTENT(INOUT) :: lx(n1,n2), lz(n1,n2), sigma(n1,n2), theta(n1,n2), tol, h
REAL,DIMENSION(:),ALLOCATABLE :: mat,mat_nnz,x,cinv_csr_val
INTEGER,DIMENSION(:),ALLOCATABLE :: irn,icn,irn_nnz,icn_nnz,cinv_csr_col,cinv_csr_ind
INTEGER::i,j,ii,jj,i1,i2,flag,lx_filesize,lz_filesize,theta_filesize,sigma_filesize,count1,nrow,nnz,operation_type,use_variable_sigma,use_rotation
REAL :: temp,SINTHETA,COSTHETA,SINTHETA2,COSTHETA2,t1,t2,lx_sc,lz_sc,sigma_sc,theta_sc,lx_norm,lz_norm,ip1,ip2,sigma_c
CHARACTER(80)::lx_file,lz_file,sigma_file,theta_file


!write(*,*)vector_in(1:10,1)
!ALLOCATE SPARSE ARRAY FOR CINV
allocate(x(n1*n2))
x(:)=0.0
allocate(irn(9*n1*n2))
allocate(icn(9*n1*n2))
allocate(mat(9*n1*n2))
irn(:)=0
icn(:)=0
mat(:)=0.0

use_rotation=1
count1=0
do i=1,n1
    do j=1,n2
        lx_sc=lx(i,j)
        lz_sc=lz(i,j)
        if (use_variable_sigma==1) then
           sigma_c=(sqrt(2*sigma(i,j)))
        else
            sigma_c=1.0
        end if
        if (use_rotation==1) then
            theta_sc=theta(i,j)
            SINTHETA=SIN(theta_sc*3.14159265359/180)
            COSTHETA=COS(theta_sc*3.14159265359/180)
            SINTHETA2=SINTHETA*SINTHETA
            COSTHETA2=COSTHETA*COSTHETA
            else
            SINTHETA=0.0
            COSTHETA=1.0
            SINTHETA2=0.0
            COSTHETA2=1.0
        end if
       lx_norm=(h/2.0)*lx_sc !(1.398*(lx_sc/h)**0.81) !100.0/lx_sc !NEED TO FIX THIS
       lz_norm=(h/2.0)*lz_sc !(1.398*(lz_sc/h)**0.81)!100.0/lz_sc !NEED TO FIX THIS
!        write(*,*)lx_sc,lz_sc,theta_sc,SINTHETA,COSTHETA,SINTHETA2,COSTHETA2,h


        ! x'=x & z'=z
        ii=i
        jj=j
        call calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
        if (flag==1) then
            count1=count1+1
            irn(count1)=i1
            icn(count1)=i2
            mat(count1)=(1/sigma_c)*(lx_norm*((1.0)/(lx_sc*h) + (2.0*lx_sc)/(h**3)) + lz_norm*((1.0)/(lz_sc*h) + (2.0*lz_sc)/(h**3)))
            !mat(count1)=(1.0/lx_sc*h)* + lx_sc*(2.0/h**3) + lz_sc*(2.0/h**3)
            vector_out(i,j)=vector_out(i,j) + vector_in(ii,jj)*mat(count1)
        end if


        ! x'=x+h & z'=z
        ii=i
        jj=j+1
        call calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
        if (flag==1) then
            count1=count1+1
            irn(count1)=i1
            icn(count1)=i2
            mat(count1)= (1/sigma_c)*(-lx_norm*(lx_sc*COSTHETA2/(h**3)) -lz_norm*(lz_sc*SINTHETA2/(h**3)))
            vector_out(i,j)=vector_out(i,j) + vector_in(ii,jj)*mat(count1)
            ! Enter value for mirror point (x'=x-h & z'=z)
            count1=count1+1
            irn(count1)=i2
            icn(count1)=i1
            mat(count1)=mat(count1-1)
            vector_out(ii,jj)=vector_out(ii,jj) + vector_in(i,j)*mat(count1)
        end if
        

        ! x'=x & z'=z+h
        ii=i+1
        jj=j
        call calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
        if (flag==1) then
            count1=count1+1
            irn(count1)=i1
            icn(count1)=i2
            mat(count1)=(1/sigma_c)*(-lx_norm*(SINTHETA2*lx_sc/(h**3)) -lz_norm*(COSTHETA2*lz_sc/(h**3)))
            vector_out(i,j)=vector_out(i,j) + vector_in(ii,jj)*mat(count1)
            ! Enter value for mirror point (x'=x & z'=z-h)
            count1=count1+1
            irn(count1)=i2
            icn(count1)=i1
            mat(count1)=mat(count1-1)
            vector_out(ii,jj)=vector_out(ii,jj) + vector_in(i,j)*mat(count1)
        end if

        ! x'=x+h & z'=z+h
        ii=i+1
        jj=j+1
        call calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
        if (flag==1) then
            count1=count1+1
            irn(count1)=i1
            icn(count1)=i2
            mat(count1)= (1/sigma_c)*(-lx_norm*((2*COSTHETA*SINTHETA*lx_sc)/(4*h**3)) + lz_norm*((2*COSTHETA*SINTHETA*lz_sc)/(4*h**3)))
            vector_out(i,j)=vector_out(i,j) + vector_in(ii,jj)*mat(count1)
           ! Enter value for mirror point (x'=x-h & z'=z-h)
            count1=count1+1
            irn(count1)=i2
            icn(count1)=i1
            mat(count1)=mat(count1-1)
            vector_out(ii,jj)=vector_out(ii,jj) + vector_in(i,j)*mat(count1)
        end if

        ! x'=x+h & z=z-h'
        ii=i+1
        jj=j-1
        call calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
        if (flag==1) then
            count1=count1+1
            irn(count1)=i1
            icn(count1)=i2
            mat(count1)=(1/sigma_c)*(lx_norm*((2*COSTHETA*SINTHETA*lx_sc)/(4*h**3)) - lz_norm*((2*COSTHETA*SINTHETA*lz_sc)/(4*h**3)))
            vector_out(i,j)=vector_out(i,j) + vector_in(ii,jj)*mat(count1)
            ! Enter value for mirror point (x'=x-h & z'=z+h)
            count1=count1+1
            irn(count1)=i2
            icn(count1)=i1
            mat(count1)=mat(count1-1)
            vector_out(ii,jj)=vector_out(ii,jj) + vector_in(i,j)*mat(count1)
        end if
        
    end do
end do


count1=0
!COUNT NUMBER OF NON-ZERO VALUES
do i=1,9*n1*n2
    if (mat(i) /= 0) then
        count1=count1+1
    end if
end do

nnz=count1
!!!write(*,*),"THE COUNT IS= ",count1
allocate(irn_nnz(nnz))
allocate(icn_nnz(nnz))
allocate(mat_nnz(nnz))

!REMOVE ZEROS VALUES FROM IRN ICN MAT
count1=0
do i=1,9*n1*n2
    if (mat(i) /= 0) then
        count1=count1+1
        irn_nnz(count1)=irn(i)
        icn_nnz(count1)=icn(i)
        mat_nnz(count1)=mat(i)
    end if
end do

call calc_n2(vector_out,n1,n2,ip1)
call calc_n2(vector_in,n1,n2,ip2)


!!!write(*,*)"IP1=",ip1," IP2=",ip2
!vector_out(:,:)=SQRT(ip2/ip1)*vector_out(:,:)

!QC SPARSE FILES
open(72,file='irn',access='direct',status='replace',recl=4*nnz)
write(72,rec=1)irn_nnz(:)
close(72)

open(72,file='icn',access='direct',status='replace',recl=4*nnz)
write(72,rec=1)icn_nnz(:)
close(72)

open(72,file='mat',access='direct',status='replace',recl=4*nnz)
write(72,rec=1)mat_nnz(:)
close(72)

if (operation_type==1) then
   allocate(cinv_csr_val(nnz))
   allocate(cinv_csr_col(nnz))
   allocate(cinv_csr_ind(nnz))
   nrow=n1*n2

   !WRITE(*,*)"ATTEMPTING TO SOLVE THE LINEAR SYSTEM TO APPLY LAPLACIAN COVARIANCE"
   !REFORMAT COO TO CSR FORMAT REQUIRED FOR CG CODE
   call coocsr (nrow,nnz,mat_nnz,irn_nnz,icn_nnz,cinv_csr_val,cinv_csr_col,cinv_csr_ind)

   call cg(cinv_csr_val,cinv_csr_ind,cinv_csr_col,nnz,nrow,nrow,vector_in,x,tol)
   !call cgmnc(cinv_csr_val,cinv_csr_ind,cinv_csr_col,nnz,nrow,nrow,vector_in,x,tol)


!!!   WRITE(*,*)"LINEAR SYSTEM SOLVED"
   call calc_n2(x,n1,n2,ip1)
   call calc_n2(vector_in,n1,n2,ip2)

!!!   write(*,*)"IP1=",ip1," IP2=",ip2
   !x(:)=SQRT(ip2/ip1)*x(:)

   call write__max(x,n1*n2)



   open(72,file='vector_out',access='direct',status='replace',recl=4*n1*n2)
   write(72,rec=1)x(:)
   close(72)

   count1=0
   do i=1,n2
      do j=1,n1
        count1=count1+1
        vector_out(j,i)=x(count1)
      end do
   end do

end if

!!!WRITE(*,*)"FINISHED USING LAPLACIAN TOOLBOX"

end subroutine lap_2d_sparse


subroutine calc2d_sparse_coords(i,j,ii,jj,n1,n2,flag,i1,i2)
implicit none

INTEGER :: i,j,ii,jj,n1,n2,flag,i1,i2

if (ii<1 .or. ii>n1) then
    flag=0
else if (jj<1 .or. jj>n2) then
    flag=0
else
    flag=1
end if

i1=(j-1)*n1+i
i2=i1+(ii-i)+(jj-j)*n1

end subroutine calc2d_sparse_coords

subroutine calc_n2(vector,n1,n2,norm_val)
implicit none

INTEGER :: i, n1,n2
REAL :: norm_val
REAL :: vector(n1*n2)
norm_val=0.0
do i=1,n1*n2
   norm_val=norm_val + vector(i)*vector(i)
end do

end subroutine calc_n2


subroutine write__max(vector,nelems)
implicit none

REAL :: vector(nelems)
REAL :: cur
INTEGER :: i,nelems
cur=0

do i=1,nelems
   if (vector(i) > cur) then
      cur=vector(i)
   end if
end do

!!!write(*,*)"MAXIMUM VALUE =",cur

end subroutine write__max




