function Gij0=Gij_PSV_r_eq_0(ks,kp,r,dr,C,vn)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj

%integration sur jj
n=2000;

d=dr/2;
% r0=logspace(-3*d,log10(d),n);
r0=linspace(-d,d,n);
i=1:n-1;
r1 =abs(r+(r0(i+1)+r0(i))/2);
% dr1=(r0(i+1)-r0(i));
g(1,:)=-vn(2)*ones(1,n-1);
g(2,:)= vn(1)*ones(1,n-1);
tmp=Gij_PSV(ks,kp,r1,g,C,n-1);
%identificacion del punto singular
%pequena trempa de evaluacion

tmp(1:2,1:2,n/2)=zeros(2,2);%r==0
%integracion
Gij0=mean(tmp,3);