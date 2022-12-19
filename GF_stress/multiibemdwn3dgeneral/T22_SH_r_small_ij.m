function t22=T22_SH_r_small_ij(coord,ii,jjpb,ki,gaussian)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj


xjj =coord.x(jjpb).';
zjj =coord.z(jjpb).';
vnx =coord.vnx(jjpb).';
vnz =coord.vnz(jjpb).';
dr  =coord.dr(jjpb).';

xi  =coord.x(ii);
zi  =coord.z(ii);
vnxi=coord.vnx(ii);
vnzi=coord.vnz(ii);

xgau=gaussian.xgau;
wgau=gaussian.wgau;

xij     = ones(gaussian.ngau,1)*xi-(xjj-(vnz.*dr)*xgau.'*0.5)*ones(1,length(xi));
zij     = ones(gaussian.ngau,1)*zi-(zjj+(vnx.*dr)*xgau.'*0.5)*ones(1,length(xi));
rij     = sqrt(xij.^2+zij.^2);
tr22    = 1i/4.*ki.* besselh(1,2,ki*rij).*(xij.*(ones(gaussian.ngau,1)*vnxi)+ zij.*(ones(gaussian.ngau,1)*vnzi))./rij;
indpb   = (rij==0);
tr22(indpb) =0;
t22     = 0.5*wgau*tr22;