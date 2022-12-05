clear all
close all
clc

nsta=40;

for k=1:nsta
   e=sprintf('out/syn_S%03d_C1.ascii',k);
   n=sprintf('out/syn_S%03d_C2.ascii',k);
   v=sprintf('out/syn_S%03d_C3.ascii',k);
   eobs2=load(e);
   nobs2=load(n);
   vobs2=load(v);

   eobs(:,1)=eobs2(:,2);%.*-1;
   nobs(:,1)=nobs2(:,2);
   vobs(:,1)=vobs2(:,2);%.*-1;

  obs = [eobs, nobs, vobs];

  e=sprintf('dat/obs_S%03d_C1',k);
  n=sprintf('dat/obs_S%03d_C2',k);
  v=sprintf('dat/obs_S%03d_C3',k);
  obsf=sprintf('dat/obs_S%03d.dat',k);

  save('-ascii',e,'eobs');
  save('-ascii',n,'nobs');
  save('-ascii',v,'vobs');
  save('-ascii',obsf,'obs');
end


