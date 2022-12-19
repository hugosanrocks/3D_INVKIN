function Gij=Gij_Pt_SVt_r_small(coord,xi,zi,jjpb,kp,gaussian,C)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xi,zi cual esta mas cercano que dr a xj,zj

xjj=coord.x(jjpb).';
zjj=coord.z(jjpb).';
vnx=coord.vnx(jjpb).';
vnz=coord.vnz(jjpb).';
dr =coord.dr(jjpb).';

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);
xij     = reshape(xi-( xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5),1,nj*gaussian.ngau);
zij     = reshape(zi-( zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5),1,nj*gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,nj*gaussian.ngau);
g       = zeros(2,nj*gaussian.ngau);
g(1,:)  = xij./rij;
g(2,:)  = zij./rij;
        
grij    = Gij_Pt_SVt(kp,rij,g,C,nj*gaussian.ngau);

% indpb   =(rij==0);
indpb   = (rij<=min(dr)/20);
grij(:,:,indpb) =Gij_Pt_SVt(kp,min(dr)/20,g(:,indpb),C,sum(indpb));

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