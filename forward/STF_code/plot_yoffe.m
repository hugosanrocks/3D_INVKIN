clear all
close all
clc

tsam=2048;

dt=25/tsam;
t=0:dt:dt*(tsam-1);
yoffe=load('yoffe.src');

for i=1:50

  plot(t,yoffe(:,i)),hold on,xlim([0,7])
  pause(0.1)
end
