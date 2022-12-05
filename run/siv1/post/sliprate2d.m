clear all
close all
clc


dt=0.25;
nsub=36*18;

%file=input('file: ')
file='../dat/modelpri.dat'
slip=load(file);
src=load('../dat/model_target.dat');
pri=load('../dat/prior_model.dat');

tsam=length(slip)/nsub;


%ti=load('dat/weight_prior.ascii');
%ti=load('weight_notmove.ascii');
%tir=reshape(ti,nsub,tsam)';
%src=load('dat/prior_model.dat');
%prinm=load('prior_yoffe.src');
%soutpri=load('source.out_e0p001_wnm');
%noreg=load('source.out_big');
%test_0p001w=load('source.out_noreg');
%test_0p001w=load('source.out_noreg');
%test_0p001=load('source.out_e0p001_1c');
%test_0p0001=load('source.out_e0p0001_1c');
%test_0p001=load('prior_2.src');

%unc=load('uncertainty.ascii');
%unc=tir.^-1;
%unc=load('uncertainty_notmove.ascii');

%timevs=load('times_vs_siv648.ascii');
%time70vs=timevs+5;
%time1p2vs=load('times_1p2vs_siv648.ascii');
%time1p2vs_r=reshape(time1p2vs,36,18)';

k=1;
for i=1:nsub
  for j=1:tsam
   s(j,i) = norm(slip(k,:));
   k = k +1;
  end
end

save('-ascii','source.out','s')


i=1;
f1=input('node to see: ')
f2=f1+1
f3=f1+2
message1=sprintf('Node # %03i',f1);
message2=sprintf('Node # %03i',f2);
message3=sprintf('Node # %03i',f3);

%source=load('source.interp_time');
%mas = zeros(17,189);
%source = [source; mas];
t=0:dt:dt*(tsam-1);

ymax=1.8;

figure(1)
subplot(311)
plot(t,src(1:tsam,f1)),hold on,
plot(t,s(:,f1),'--r'),hold on,
plot(t,pri(1:tsam,f1),'*-k')
xlabel('Time (sec)'),ylabel('Slip rate (m/s)')
title(message1)
legend('SIV1','Inv','Prior');
set(gca,'Fontsize',15);
ylim([0,ymax])

figure(1)
subplot(312)
plot(t,src(1:tsam,f2)),hold on,
plot(t,s(:,f2),'--r')
plot(t,pri(1:tsam,f2),'*-k')
xlabel('Time (sec)'),ylabel('Slip rate (m/s)')
title(message2)
%legend('Target','Inversion');
set(gca,'Fontsize',15);
ylim([0,ymax])

figure(1)
subplot(313)
plot(t,src(1:tsam,f3)),hold on,
plot(t,s(:,f3),'--r')
plot(t,pri(1:tsam,f3),'*-k')
title(message3)
xlabel('Time (sec)'),ylabel('Slip rate (m/s)')
%legend('Target','Inversion');
set(gca,'Fontsize',15);
ylim([0,ymax])


print('../graphics/sliprate_comparison.pdf','-dpdf')


