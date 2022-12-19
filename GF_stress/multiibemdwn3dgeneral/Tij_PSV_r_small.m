function Tij=Tij_PSV_r_small(coord,ii,jjpb,ks,kp,gaussian,C)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj


xjj     = coord.x(jjpb).';
zjj     = coord.z(jjpb).';
vnx     = coord.vnx(jjpb).';
vnz     = coord.vnz(jjpb).';
dr      = coord.dr(jjpb).';

xi      = coord.x(ii);
zi      = coord.z(ii);
vn(1,1) = coord.vnx(ii);
vn(2,1) = coord.vnz(ii);

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);

xij     = reshape(xi-(xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5),1,nj*gaussian.ngau);
zij     = reshape(zi-(zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5),1,nj*gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,nj*gaussian.ngau);
g       = zeros(2,nj*gaussian.ngau);

g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

n=nj;

Trij=Tij_PSV(ks,kp,rij,g,C,vn);
indpb   =(rij==0);
Trij(:,:,indpb) =0;
Trij=reshape(Trij,2,2,nj,gaussian.ngau);


Tij=zeros(2,2,n);
if nj~=1
    for i=1:2
        for j=1:2
            Tij(i,j,:)	= sum(0.5*(ones(nj,1)*wgau).*squeeze(Trij(i,j,:,:)),2);
        end
    end
else
    for i=1:2
        for j=1:2
            Tij(i,j,:)	= sum(0.5*wgau.'.*squeeze(Trij(i,j,:,:)));
        end
    end
end