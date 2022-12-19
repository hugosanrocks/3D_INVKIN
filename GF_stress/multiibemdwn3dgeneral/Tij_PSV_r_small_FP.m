function Tij=Tij_PSV_r_small_FP(coord,xs,zs,jjpb,ks,kp,gaussian,C)
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
g       = zeros(2,nj*gaussian.ngau);

g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

n=nj;

Trij=Tij_PSV(ks,kp,rij,g,C,vn);
indpb   =(rij==0);
Trij(:,:,indpb) =0;
if max(indpb)==0
    [~,indmr]=sort(rij);
    if max(abs(Trij(:,:,indmr(1))))>3*max(abs(Trij(:,:,indmr(2))))
        Trij(:,:,indmr(1)) = Trij(:,:,indmr(2))*2.*sign(Trij(:,:,indmr(1))).*sign(Trij(:,:,indmr(2)));
    end
%     tmp=diff(sign(real(Trij(1,2,indmr))));
%     indok=find(tmp~=0,1,'last')+1;
%     Trij(1,2,indmr(1:indok)) = 0;
%     Trij(2,1,indmr(1:indok)) = 0;
%     if indmr(1)>1 && indmr(1)<nj
%         Trij(:,:,indmr(1)) = .5*(Trij(:,:,indmr(1)+1)+Trij(:,:,indmr(1)-1));
%     else
%         Trij(1,2,indmr(1:indok)) = 0;
%         Trij(2,1,indmr(1:indok)) = 0;
%     end
end
Trij=reshape(Trij,2,2,nj,gaussian.ngau);
for i=1:gaussian.ngau
    tmp=.5*mean(abs(Trij(1,2,:,i)));
    indpb=find(abs(Trij(1,2,:,i))>tmp);
    Trij(1,2,indpb,i)=0;
    Trij(2,1,indpb,i)=0;
end

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