program testLoad

IMPLICIT NONE

   INTEGER, PARAMETER :: dp = KIND(1D0)
   INTEGER :: funit, io_stat
   COMPLEX(dp), allocatable :: s(:,:)
   real(dp), allocatable :: rawMat(:,:)
   integer :: i,j,k,f,nfrec,nrow,ncol,nFault,nSta
   COMPLEX(dp), parameter :: UI = cmplx(0.0,1.0,dp)

   print*, 'leer información del tamaño de la variable'
   OPEN(123,FILE = "info.txt",FORM="FORMATTED")
   read(123,*) !Number of frequencies:
   read(123,*) nfrec
   nfrec = 5521
   read(123,*) !Number of rows:
   read(123,*) nrow
   read(123,*) !Number of columns:
   read(123,*) ncol
   read(123,*) !Number of Fault points(receibers):
   read(123,*) nFault
   read(123,*) !Number of Stations (sources):
   read(123,*) nSta
   close(123)
   print*,nfrec,nrow,ncol,nFault,nSta

   ! print True o False si la cantidad de registros es
   ! la cantidad de puntos de en la falla * estaciones
   print*, (nFault*nSta == nrow) !debe ser T

   ! dimensionar variables
   allocate(rawMat(nrow,ncol)) ! variable auxiliar
   allocate(s(nrow,nfrec))     ! el vector con los esfuerzos

   OPEN(NEWUNIT = funit, FILE = 'SIGMA_XX_x', STATUS = "OLD", ACCESS = "STREAM", FORM = "UNFORMATTED", IOSTAT = io_stat)
   READ(funit, IOSTAT = io_stat) rawMat
   CLOSE(funit)

   ! en rawMat la parte real e imaginaria están como matrices adjacentes
   ! las copiamos a 's' como una variable compleja:
   !do i=1,nrow
   !  s(i,:) = rawMat(i,1:nfrec) + UI * rawMat(i,nfrec+1:ncol)
   !end do

   ! fin

! ver los valores para todos los registros en una frecuencia dada
   f = 5 ! <-- indice de la frecuencia mostrada
   i=1;
   !do k=1,nFault
   !do j=1,nSta
   !write(6,'(a,I0,a,I0,a,I0,a,EN12.3,1x,EN12.3,a)') &
   !     'registro',i,'_sta',j,'fault',k,'=   ',real(s(i,f)),imag(s(i,f)),'i'
   !i=i+1
   !end do
   !end do
   do i=1,ncol
   write(33,*) rawMat(:,i)
   enddo
end program testLoad
