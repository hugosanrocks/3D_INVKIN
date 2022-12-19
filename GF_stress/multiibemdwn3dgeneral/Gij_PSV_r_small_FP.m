function Gij=Gij_PSV_r_small_FP(coord,xs,zs,j,ks,kp,gaussian,C)
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
g       = zeros(2,nj*gaussian.ngau);
g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

grij    = Gij_PSV(ks,kp,rij,g,C,nj*gaussian.ngau);

indpb   =(rij==0);
grij(:,:,indpb) =0;

if max(indpb)==0
    [~,indmr]=sort(rij);
    if max(abs(grij(:,:,indmr(1))))>3*max(abs(grij(:,:,indmr(2))))
        grij(:,:,indmr(1)) = grij(:,:,indmr(2))*2.*sign(grij(:,:,indmr(1))).*sign(grij(:,:,indmr(2)));
    end
%     tmp=diff(sign(real(grij(1,2,indmr))));
%     indok=find(tmp~=0,1,'last');
%     grij(1,2,indmr(1:indok)) = 0;
%     grij(2,1,indmr(1:indok)) = 0;
end

grij    = reshape(grij,2,2,nj,gaussian.ngau);
Gij     = zeros(2,2,nj);

if nj~=1
    Gij(1,1,:)	= sum(0.5*(ones(nj,1)*wgau).*squeeze(grij(1,1,:,:)),2);
    Gij(2,2,:)	= sum(0.5*(ones(nj,1)*wgau).*squeeze(grij(2,2,:,:)),2);
    Gij(1,2,:)	= sum(0.5*(ones(nj,1)*wgau).*squeeze(grij(1,2,:,:)),2);
    Gij(2,1,:)	= Gij(1,2,:);
else
    Gij(1,1,:)	= sum(0.5*wgau.'.*squeeze(grij(1,1,:,:)));
    Gij(2,2,:)	= sum(0.5*wgau.'.*squeeze(grij(2,2,:,:)));
    Gij(1,2,:)	= sum(0.5*wgau.'.*squeeze(grij(1,2,:,:)));
    Gij(2,1,:)	= Gij(1,2,:);
end
