function [S1,S2]=S_PSV_r_small_2(coordf,xr,zr,ks,kp,C,fij,gaussian)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xr,zr cual esta mas cercano que dr a xj,zj

xjj     = coordf.x;
zjj     = coordf.z;
vnx     = coordf.vnx;
vnz     = coordf.vnz;
dr      = coordf.dr;

xgau    = gaussian.xgau;
wgau    = gaussian.wgau;

xij     = reshape(xr-(xjj*ones(1,gaussian.ngau)-(vnz.*dr)*xgau*0.5),1,gaussian.ngau);
zij     = reshape(zr-(zjj*ones(1,gaussian.ngau)+(vnx.*dr)*xgau*0.5),1,gaussian.ngau);
rij     = reshape(sqrt(xij.^2+zij.^2),1,gaussian.ngau);
g       = zeros(2,gaussian.ngau);

g(1,:)  = xij./rij;
g(2,:)  = zij./rij;

[S1,S2] = S_PSV(rij,g,ks,kp,C,fij);

indpb       = (rij==0);

S           = S1;
S(:,:,indpb)= 0;
S           = reshape(S,2,2,gaussian.ngau);
S1          = zeros(2,2);
i=1;
for j=1:2
    S1(i,j)	= sum(0.5*wgau.'.*squeeze(S(i,j,:)));
end
j=2;i=2;
S1(i,j)     = sum(0.5*wgau.'.*squeeze(S(i,j,:)));

S           = S2;
S(:,:,indpb)= 0;
S           = reshape(S,2,2,gaussian.ngau);
S2          = zeros(2,2);
i=1;
for j=1:2
    S2(i,j)	= sum(0.5*wgau.'.*squeeze(S(i,j,:)));
end
j=2;i=2;
S2(i,j)     = sum(0.5*wgau.'.*squeeze(S(i,j,:)));
