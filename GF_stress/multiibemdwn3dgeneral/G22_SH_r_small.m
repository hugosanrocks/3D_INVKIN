function gn22=G22_SH_r_small(coord,xi,zi,jjpb,ki,gaussian,mu)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj


xjj=coord.x(jjpb).';
zjj=coord.z(jjpb).';
vnx=coord.vnx(jjpb).';
vnz=coord.vnz(jjpb).';
dr =coord.dr(jjpb).';

% vtx= vnz; %vecteur tangent
% vtz=-vnx;
% 
% xjj0=xjj+vtx*dr/2;
% zjj0=zjj+vtz*dr/2;
% 
% xjj1=xjj-vtx*dr/2;
% zjj1=zjj-vtz*dr/2;
% 
% n=400;
% i=1:n-1;
% 
% xint=linspace(xjj0,xjj1,n);
% xint1=0.5*(xint(i)+xint(i+1));
% 
% zint=linspace(zjj0,zjj1,n);
% zint1=0.5*(zint(i)+zint(i+1));
% 
% xij0 = xi-xint1;
% zij0 = zi-zint1;
% rij0= sqrt(xij0.^2+zij0.^2);
% drij= sqrt((xint(i)-xint(i+1)).^2+(zint(i)-zint(i+1)).^2);
% gn220=mean(G22_SH(ki,rij0,mu));

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);
xij     = xi-( xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5);
zij     = zi-( zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5);
rij     = sqrt(xij.^2+zij.^2);
gr22    = G22_SH(ki,rij,mu);
indpb   =(rij==0);
gr22(indpb) =0;
gn22    = sum(0.5*(ones(nj,1)*wgau).*gr22,2);