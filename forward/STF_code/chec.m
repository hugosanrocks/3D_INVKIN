clear all
close all
clc

dt=25/2048;
t=0:dt:dt*2046;
s=load('fort.33');
plot(t,s),xlim([0,4])
