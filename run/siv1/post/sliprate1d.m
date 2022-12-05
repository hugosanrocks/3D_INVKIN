clear all
close all
clc

nsub=36*18;
dt=0.25;
filein='../dat/modelpri1d.dat'

%Load the file and detect time samples
slip=load(filein);
tsam=length(slip)/nsub

%load prior model
fileprior='../dat/prior_model.dat'
prior=load(fileprior);

%load target of inversion
filesource='../dat/model_target.dat'
target=load(filesource);

%build time axis to plot
tim=tsam*dt-dt
t=0:dt:dt*(tsam-1);

%reorganize the solution
k=1;
for i=1:nsub
  for j=1:tsam
   sliprate(j,i) = slip(k);
   k=k+1;
  end
end

%save inversion results
fileout='source.out'
save('-ascii',fileout,'sliprate')

%IF YOU WANT TO PLOT SOME NODE HISTORY

%for i=1:36
i=1;
f1=input('node to see: ')
f2=f1+1
f3=f1+2
message1=sprintf('Node # %03i',f1);
message2=sprintf('Node # %03i',f2);
message3=sprintf('Node # %03i',f3);

%ma = max(max(p1));
ma=1.8;
tfin=9;

figure(i)
subplot(3,1,1)
%target
a=plot(t,target(1:tsam,f1),'k'),hold on,
%inversion
b=plot(t,sliprate(1:tsam,f1),'.-r')
%prior model
c=plot(t,prior(1:tsam,f1),'.-b')
%other plots
%d=plot(t,tri(1:tsam,f1),'*-b')
%e=plot(t,sourcep(1:tsam,f1),'+-b')
%plot(t,pre(:,f1),'--y');
legend([a,b,c],'SIV1','Inv','Prior')
ylim([0,ma])
xlim([0,tfin])
title(message1)

subplot(3,1,2)
%target
plot(t,target(1:tsam,f2),'k'),hold on,
%inversion
plot(t,sliprate(:,f2),'.-r')
hold on
plot(t,prior(1:tsam,f2),'.-b')
%other plots
%plot(t,tri(1:tsam,f2),'*-b')
%plot(t,sourcep(1:tsam,f2),'+-b')
%plot(t,pre(:,f2),'--y');
%plot(t,p1(:,f2),'*-b'),plot(t,p2(:,f2),'*-b')
ylim([0,ma])
xlim([0,tfin])
title(message2)

subplot(3,1,3)
%target
plot(t,target(1:tsam,f3),'k'),hold on,
%inversion
plot(t,sliprate(:,f3),'.-r')
hold on
plot(t,prior(1:tsam,f3),'.-b')
%plot(t,tri(1:tsam,f3),'*-b')
%plot(t,sourcep(1:tsam,f3),'+-b')
%plot(t,pre(:,f3),'--y');
%plot(t,p1(:,f3),'*-b'),plot(t,p2(:,f3),'*-b')
ylim([0,ma])
xlim([0,tfin])
title(message3)
print('../graphics/sliprate_comparison.pdf','-dpdf')

%end

