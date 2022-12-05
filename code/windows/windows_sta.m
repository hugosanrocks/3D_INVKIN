clear all
close all
clc


w1 = 5:2:39;
w2 = 11:5:96;

dt=0.1;
nsamp=351;
time=input('time: ');
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


