function [S_fx,S_fy,S_fz]=S_3D_r_small(coordf,jjpb,xr,yr,zr,ks,kp,fij,gaussian)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xr,zr cual esta mas cercano que dr a xj,zj

xfc0    = coordf.x0(jjpb).';
thfc    = coordf.th(jjpb).';
nbptc   = coordf.nbptc(jjpb).';

zfc     = coordf.z(jjpb).';

vnx     = coordf.vnx(jjpb).';
vnz     = coordf.vnz(jjpb).';
vnx0    = vnx./cos(thfc);
ind     = (mod(thfc-pi/2,pi)==0);
vnx0(ind)= coordf.vny(jjpb(ind)).'./sin(thfc(ind));
    
drxz    = coordf.drxz(jjpb).';

ngau    = gaussian.ngau;
xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(xfc0);

%construction des coordonnees des points de Gauss
xij     = zeros(nj,ngau^2);
yij     = zeros(nj,ngau^2);
zij     = zeros(nj,ngau^2);
wij     = zeros(nj,ngau^2);
k       = 0;
% attention il faut rajouter la correction des aires

% figure;hold on
% plot3(coordf.x,coordf.y,coordf.z,'.')

for i=1:ngau
    %boucle sur le vertical, plan xz
    % pour le calcul des rayons initiaux
    xfi     = xfc0-vnz .*drxz*xgau(i)*0.5;
    zfi     = zfc +vnx0.*drxz*xgau(i)*0.5;
    for j=1:ngau
        k=k+1;
        %boucle sur le rayon ds le plan xy
        % on recalcule les positions des points de gauss 
        % et leur normales(car elements courbes)
        % on dispose les points de Gauss sur l element courbe comme si
        % deplies (le long de r)
        dthgauss= 2*pi./nbptc*xgau(j)*.5;
        xij(:,k)= xr-xfi.*cos(thfc+dthgauss);
        yij(:,k)= yr-xfi.*sin(thfc+dthgauss);
        zij(:,k)= zr-zfi;
        wij(:,k)= wgau(i)*wgau(j)*xfi./xfc0;
%         plot3(xfi(1).*cos(thfc(1)+dthgauss(1)),xfi(1).*sin(thfc(1)+dthgauss(1)),zfi(1),'r.')
    end
end
xij     = reshape(xij,1,nj*ngau^2);
yij     = reshape(yij,1,nj*ngau^2);
zij     = reshape(zij,1,nj*ngau^2);

rij     = sqrt(xij.^2+yij.^2+zij.^2);

g       = zeros(3,nj*ngau^2);
g(1,:)  = xij./rij;
g(2,:)  = yij./rij;
g(3,:)  = zij./rij;




[S_fx0,S_fy0,S_fz0] = S_3D(rij,g,ks,kp,fij);

indpb               = (rij==0);
S_fx0(:,:,indpb)    = 0;
S_fx0               = reshape(S_fx0,3,3,nj,ngau^2);
S_fy0(:,:,indpb)    = 0;
S_fy0               = reshape(S_fy0,3,3,nj,ngau^2);
S_fz0(:,:,indpb)    = 0;
S_fz0               = reshape(S_fz0,3,3,nj,ngau^2);


S_fx=zeros(3,3,nj);
if nj~=1
    for i=1:3
        for j=1:3
            S_fx(i,j,:)	= sum(0.25*wij.*squeeze(S_fx0(i,j,:,:)),2);
        end
    end
else
    for i=1:3
        for j=1:3
            S_fx(i,j,:)	= sum(0.25*wij.'.*squeeze(S_fx0(i,j,:,:)));
        end
    end
end

S_fy=zeros(3,3,nj);
if nj~=1
    for i=1:3
        for j=1:3
            S_fy(i,j,:)	= sum(0.25*wij.*squeeze(S_fy0(i,j,:,:)),2);
        end
    end
else
    for i=1:3
        for j=1:3
            S_fy(i,j,:)	= sum(0.25*wij.'.*squeeze(S_fy0(i,j,:,:)));
        end
    end
end

S_fz=zeros(3,3,nj);
if nj~=1
    for i=1:3
        for j=1:3
            S_fz(i,j,:)	= sum(0.25*wij.*squeeze(S_fz0(i,j,:,:)),2);
        end
    end
else
    for i=1:3
        for j=1:3
            S_fz(i,j,:)	= sum(0.25*wij.'.*squeeze(S_fz0(i,j,:,:)));
        end
    end
end