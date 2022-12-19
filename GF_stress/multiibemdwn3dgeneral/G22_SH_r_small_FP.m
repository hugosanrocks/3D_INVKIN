function gn22=G22_SH_r_small_FP(coord,xs,zs,j,ki,gaussian,mu)
% trata el calculo de la contribution de una fuente en xs,zs
% sobre un segmento de un punto de colocacion xj,zj cual esta mas cercano
% que dr a xs,zs

xjj=coord.x(j).';
zjj=coord.z(j).';
vnx=coord.vnx(j).';
vnz=coord.vnz(j).';
dr =coord.dr(j).';

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);
xij     = reshape(( xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5)-xs,1,nj*gaussian.ngau);
zij     = reshape(( zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5)-zs,1,nj*gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,nj*gaussian.ngau);
        
gr22    = G22_SH(ki,rij,mu);
gr22    =reshape(gr22,nj,gaussian.ngau);

indpb   =(rij==0);
gr22(indpb) =0;
gn22    = sum(0.5*(ones(nj,1)*wgau).*gr22,2);
