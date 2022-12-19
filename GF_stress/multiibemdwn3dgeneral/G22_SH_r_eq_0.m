function gn22=G22_SH_r_eq_0(ki,dr,mu)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj
% igual que antiguo Greenex
% %integration sur jj
% n=2000;
% 
% d=dr/2;
% % r0=logspace(-3*d,log10(d),n);
% r0=linspace(-d,d,n);
% i=1:n-1;
% r1 =abs((r0(i+1)+r0(i))/2);
% % dr1=(r0(i+1)-r0(i));
% indtmp=find(abs(r1)<1e-10);
% r1(indtmp)=r1(indtmp+1);
% gn22=mean(G22_SH(ki,r1,mu));

ks  = ki;
geu	= 0.57721566;
i2pi= 1i*2/pi;
fac1= 1-i2pi*(geu-1);
fac2= 1-i2pi*(geu-1-1/3);
logs= log(ks*dr/4);
h0S = (fac1-i2pi*logs) - dr.^2*ks.^2/48*(fac2-i2pi*logs);
gn22= 1/(4*1i*mu).*h0S;
