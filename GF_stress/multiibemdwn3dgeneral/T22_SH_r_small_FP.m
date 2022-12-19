function t22=T22_SH_r_small_FP(coord,xs,zs,jjpb,ks,gaussian)
% trata del calculo de la contribution de una fuente en xs,zs
% sobre un segmento de un punto de colocacion xj,zj cual esta mas cercano
% que dr a xs,zs

xjj     = coord.x(jjpb).';
zjj     = coord.z(jjpb).';
vnx     = coord.vnx(jjpb).';
vnz     = coord.vnz(jjpb).';
dr      = coord.dr(jjpb).';

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);

vn      = zeros(2,nj*gaussian.ngau);
vn(1,:) = reshape(vnx*ones(1,gaussian.ngau),1,nj*gaussian.ngau);
vn(2,:) = reshape(vnz*ones(1,gaussian.ngau),1,nj*gaussian.ngau);

xij     = reshape((xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5)-xs,1,nj*gaussian.ngau);
zij     = reshape((zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5)-zs,1,nj*gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,nj*gaussian.ngau);
tr22    = 1i/4.*ks.* besselh(1,2,ks*rij).*(xij.*vn(1,:)+ zij.*vn(2,:))./rij;
indpb   = (rij==0);
tr22(indpb) =0;
tr22=reshape(tr22,nj,gaussian.ngau);
t22     = sum(0.5*(ones(nj,1)*wgau).*tr22,2);
