clear all
close all
clc


w1 = 5:2:39;
w2 = 11:5:96;

dt=0.25;
nsamp=19;
nsta = 40;
ncomp = 3;
water = 0.01;

stacomp = nsta*ncomp;

%w=load('dat/weights.dat');
%w(:,:) =1;
for k=1:nsta
   eo=sprintf('dat/obs_S%03d_C1',k);
   no=sprintf('dat/obs_S%03d_C2',k);
   vo=sprintf('dat/obs_S%03d_C3',k);
   es=sprintf('out/syn_S%03d_C1.ascii',k);
   ns=sprintf('out/syn_S%03d_C2.ascii',k);
   vs=sprintf('out/syn_S%03d_C3.ascii',k);
   eobs(:,k)=load(eo);
   nobs(:,k)=load(no);
   vobs(:,k)=load(vo);
   esyn2c=load(es);
   nsyn2c=load(ns);
   vsyn2c=load(vs);
end

snaps=[3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39];
n=1:nsta;

for k=1:nsamp
   time=(snaps(k)*dt)-dt;
   message2=sprintf('window_%02.2f',time)
   window(:,k+2)=load(message2);
end
window(:,1)=n';
window(:,2)=0;

for i=1:nsta
message3=sprintf('%02i %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f %02.2f',window(i,:))
fileout=fopen('windows.dat','a');
fprintf(fileout,'%s\n',message3);
end

%save('-ascii','windows.dat','window');
