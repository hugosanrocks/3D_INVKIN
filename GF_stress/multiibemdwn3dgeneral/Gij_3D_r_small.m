function Gij=Gij_3D_r_small(coord,xr,yr,zr,jjpb,ks,kp,gaussian,C)
% trata el calculo de la contribution de un elemento centrado en xj,yj,zj
% sobre un punto de colocacion xi,yi,zi cual esta mas cercano que dr a xj,yj,zj

xfc0    = coord.x0(jjpb).';
thfc    = coord.th(jjpb).';
nbptc   = coord.nbptc(jjpb).';

zfc     = coord.z(jjpb).';

vnx     = coord.vnx(jjpb).';
vnz     = coord.vnz(jjpb).';
vnx0    = vnx./cos(thfc);
ind     = (mod(thfc-pi/2,pi)==0);
vnx0(ind)= coord.vny(jjpb(ind)).'./sin(thfc(ind));
    
drxz    = coord.drxz(jjpb).';

ngau    = gaussian.ngau;
xgau    = gaussian.xgau;
wgau    = gaussian.wgau;
nj      = length(jjpb);

%construction des coordonnees des points de Gauss
xij = zeros(nj,ngau^2);
yij = zeros(nj,ngau^2);
zij = zeros(nj,ngau^2);
wij = zeros(nj,ngau^2);
k   = 0;
% attention il faut rajouter la correction des aires

% figure;hold on
% plot3(coord.x,coord.y,coord.z,'.')

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

grij    = Gij_3D(ks,kp,rij,g,C,nj*ngau^2);
indpb   = (rij==0);
grij(:,:,indpb) =0;   % <- suspicious TODO
grij    = reshape(grij,3,3,nj,ngau^2);


Gij     = zeros(3,3,nj);
if nj~=1
    for i=1:3
        for j=i:3
            Gij(i,j,:)	= sum(0.25*wij.*squeeze(grij(i,j,:,:)),2);
            Gij(j,i,:)	= Gij(i,j,:);
        end
    end
else
    for i=1:3
        for j=i:3
            Gij(i,j,:)	= sum(0.25*wij.'.*squeeze(grij(i,j,:,:)));
            Gij(j,i,:)	= Gij(i,j,:);
        end
    end
end