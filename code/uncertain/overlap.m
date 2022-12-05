clear all
close all
clc

tsam=60;
dt=0.25;
t=0:dt:dt*(tsam-1);
tw=2;
sam_w=floor(tw/dt)+1;

ini=floor(sam_w/2);
fin=tsam-ini;

