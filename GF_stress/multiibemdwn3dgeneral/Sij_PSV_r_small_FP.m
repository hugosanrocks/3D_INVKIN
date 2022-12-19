function Sij = Sij_PSV_r_small_FP(coord,xs,zs,j,ks,kp,gaussian,C)
% Función de Green de esfuerzos en el campo cercano
% trata el calculo de la contribution de una fuente en xs,zs
% sobre un segmento de un punto de colocacion xj,zj cual esta mas cercano
% que dr a xs,zs

xjj     = coord.x(j).';
zjj     = coord.z(j).';
vnx     = coord.vnx(j).';
vnz     = coord.vnz(j).';
dr      = coord.dr(j).';

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xjj);
xij     = reshape(( xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5)-xs,1,nj*gaussian.ngau);
zij     = reshape(( zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5)-zs,1,nj*gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,nj*gaussian.ngau);
g       = zeros(2,nj*gaussian.ngau);
g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

srij    = Sij_PSV(rij,g,ks,kp,C); % (sxx|sxz|szz,Fza x|z,:)

indpb   = (rij==0);
srij(:,:,indpb) =0;

if max(indpb)==0
    [~,indmr]=sort(rij);
    if max(abs(srij(:,:,indmr(1))))>3*max(abs(srij(:,:,indmr(2))))
        srij(:,:,indmr(1)) = srij(:,:,indmr(2))*2.*sign(srij(:,:,indmr(1))).*sign(srij(:,:,indmr(2)));
    end
end

srij    = reshape(srij,3,2,nj,gaussian.ngau);
Sij     = zeros(3,2,nj);

if nj~= 1
  Sij(1,1,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(1,1,:,:)),2);
  Sij(2,1,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(2,1,:,:)),2);
  Sij(3,1,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(3,1,:,:)),2);
  Sij(1,2,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(1,2,:,:)),2);
  Sij(2,2,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(2,2,:,:)),2);
  Sij(3,2,:) = sum(0.5*(ones(nj,1)*wgau).*squeeze(srij(3,2,:,:)),2);
else
  Sij(1,1,:) = sum(0.5*wgau.'.*squeeze(srij(1,1,:,:)));
  Sij(2,1,:) = sum(0.5*wgau.'.*squeeze(srij(2,1,:,:)));
  Sij(3,1,:) = sum(0.5*wgau.'.*squeeze(srij(3,1,:,:)));
  Sij(1,2,:) = sum(0.5*wgau.'.*squeeze(srij(1,2,:,:)));
  Sij(2,2,:) = sum(0.5*wgau.'.*squeeze(srij(2,2,:,:)));
  Sij(3,2,:) = sum(0.5*wgau.'.*squeeze(srij(3,2,:,:)));
end
end