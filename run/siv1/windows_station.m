clear all
close all
clc

%time sampling of seismograms
dt=0.1;
%number of samples
nsamp=351;
%time window limit for the source
time=input('time: ');
%number of stations
nsta = 40;
%components
ncomp = 3;
%percentage of difference allowed
%between synthetic and observations
water = 0.01;

stacomp = nsta*ncomp;

%read observations and synthetics
%and limit the time window where
%their difference is meaningful
for k=1:nsta
   fileobs=sprintf('dat/obs_S%03d.dat',k);
   es=sprintf('out/syn_S%03d_C1.ascii',k);
   ns=sprintf('out/syn_S%03d_C2.ascii',k);
   vs=sprintf('out/syn_S%03d_C3.ascii',k);
   obs=load(fileobs);
   eobs(:,k)=obs(:,1);
   nobs(:,k)=obs(:,2);
   vobs(:,k)=obs(:,3);
   esyn2c=load(es);
   nsyn2c=load(ns);
   vsyn2c=load(vs);
   esyn(:,k)=esyn2c(:,2);
   nsyn(:,k)=nsyn2c(:,2);
   vsyn(:,k)=vsyn2c(:,2);
   me(k)=(max(abs(eobs(:,k)))*water)^2;
   mn(k)=(max(abs(nobs(:,k)))*water)^2;
   mv(k)=(max(abs(vobs(:,k)))*water)^2;

   eres(:,k) = eobs(:,k) - esyn(:,k);
   nres(:,k) = nobs(:,k) - nsyn(:,k);
   vres(:,k) = vobs(:,k) - vsyn(:,k);
   eres2(:,k)=eres(:,k).^2;
   nres2(:,k)=nres(:,k).^2;
   vres2(:,k)=vres(:,k).^2;
  
   ie = 1;
   in = 1;
   iv = 1; 

   while ( (eres2(ie,k) < me(k)) && (ie < nsamp))
     ie = ie+1;
   end
   inde(k) = ie;

   while ( (nres2(in,k) < mn(k)) && (in < nsamp))
     in = in+1;
   end
   indn(k) = in;

   while ( (vres2(iv,k) < mv(k)) && (iv < nsamp))
     iv = iv+1;
   end
   indv(k) = iv;


end

sta=1;%input('sta: ');

y=[-1,1];
x=[inde(sta),inde(sta)];
plot(eobs(:,sta),'k'),hold on
plot(esyn(:,sta),'.-r')
plot(eres2(:,sta).*50,'b','linewidth',0.2)
plot(x,y,'g')
ylim([min(eobs(:,sta)),max(eobs(:,sta))])

indet = (inde.*dt) - dt;
indnt = (indn.*dt) - dt;
indvt = (indv.*dt) - dt;
indet = indet';
indnt = indnt';
indvt = indvt';
for i=1:nsta
 minind(i,1) = min([indet(i),indnt(i),indvt(i)]);
end

y2=[minind(sta) minind(sta)];

message=sprintf('window_%02.2f',time)
save('-ascii',message,'minind')


