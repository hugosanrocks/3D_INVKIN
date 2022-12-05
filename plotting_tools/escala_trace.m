clear all
close all
clc

i=input('scale for station:');

    file=sprintf('out/obs_S%03d_C1.a',i);
    syn=load(file);
    maxa=max(abs(syn(:,2)));
    tmax=max(syn(:,1));
    esca1 = [tmax-3 maxa*0.5; tmax-1 maxa*0.5];
    esca2 = [tmax-3 maxa*0.25; tmax-1 maxa*0.25];
    save('-ascii','graphics/escala1.dat','esca1')
    save('-ascii','graphics/escala2.dat','esca2')

