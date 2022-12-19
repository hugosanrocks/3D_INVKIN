function t22=T22_SH_r_small(coord,ii,jjpb,ki,gaussian)
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
nj  =length(xjj);

xij     = xi-(xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5);
zij     = zi-(zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5);
rij     = sqrt(xij.^2+zij.^2);
tr22    = 1i/4.*ki.* besselh(1,2,ki*rij).*(xij.*vnxi+ zij.*vnzi)./rij;
indpb   = (rij==0);
tr22(indpb) =0;
t22     = sum(0.5*(ones(nj,1)*wgau).*tr22,2);