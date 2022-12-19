function [Sxy,Szy]=Sxzy_SH_r_small(coord,xr,zr,ki,gaussian)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xr,zr cual esta mas cercano que dr a xj,zj


xjj =coord.x;
zjj =coord.z;
vnx =coord.vnx;
vnz =coord.vnz;
dr  =coord.dr;

xgau=gaussian.xgau;
wgau=gaussian.wgau;

xij     = xr-(xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5);
zij     = zr-(zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5);
rij     = sqrt(xij.^2+zij.^2);

tmp    = 1i/4.*ki.* besselh(1,2,ki*rij);
Sxy    = tmp.*xij./rij;
Szy    = tmp.*zij./rij;
indpb  = (rij<min(dr)*2);

Sxy(indpb) =0;
Szy(indpb) =0;

nj      = length(xjj);
Sxy     = sum(0.5*(ones(nj,1)*wgau).*Sxy,2);
Szy     = sum(0.5*(ones(nj,1)*wgau).*Szy,2);
