clear all
close all
clc

%time sampling
dt=0.25;
%number of windows
nsamp=9;
%stations
nsta = 40;
%components
ncomp = 3;

system('rm windows.dat')

stacomp = nsta*ncomp;


snaps=[3 5 7 9 13 19 27 35 39];
n=1:nsta;

for k=1:nsamp
   time=(snaps(k)*dt)-dt;
   message2=sprintf('window_%02.2f',time)
   window(:,k+2)=load(message2);
end
window(:,1)=n';
window(:,2)=0;

%write time limits for each time window
for i=1:nsta
message3=sprintf('%02i %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f',window(i,:));
fileout=fopen('windows.dat','a');
fprintf(fileout,'%s\n',message3);
end

%save file in correct directory
system('mv windows.dat dat/windows.info');
display('Time-window limits on dat/windows.info')

