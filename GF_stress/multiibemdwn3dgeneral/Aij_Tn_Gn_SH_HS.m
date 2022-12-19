function [t22,gn22]=Aij_Tn_Gn_SH_HS(j,coord,para)
% funcion calculando traciones normales y desplazamientos
% es una funcion similar a Aij_Tn_Gn_SH pero para un semi-espacio(Half-Space)
% considerando fuentes virtuales imagen
% Asi da los coeficientes de la matriz A que corresponden a la continuidad
% de los esfuerzos normales 
% j : indice del punto de la fuente virtual ponctual que tiene una
% contribucion sobre el conjunto de los elementos de la superficie i

%%%%%%%%%%%%%
% unpacking %
%%%%%%%%%%%%%
indpt           = 1:coord.nbpt;
indi            = coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1
indir           = indpt(indi);%indice de estos puntos

%coordenadas de los elementos receptores
xr              = coord.x(indi);
zr              = coord.z(indi);
vnx             = coord.vnx(indi);
vnz             = coord.vnz(indi);
dr              = coord.dr(indi);
% jphi            = coord.phi(indi,1);
drj             = coord.dr(j);
% t22             = zeros(1,coord.nbeq);
% t221            = zeros(1,coord.nbeq);
% gn22            = zeros(1,coord.nbeq);
% gn221           = zeros(1,coord.nbeq);

gaussian           = para.gaussian;

% %indice (logic y natural) de los puntos de colocacion perteneciendo a m
% jjx             = coord.indm(m).ind;
% jr              = j(jjx);
% 
% %indice (logic y natural) de los phi que hay que contemplar (columnas
% %de la matriz)
% jjphi           = false(coord.nbeq,1);
% jjphi(phi(j(jjx),m)) = true(1);
% jj              = jphi(jjphi);

ksi             = para.reg(1).sub(1).ksi;
mu              = para.reg(1).sub(1).Ci(6,6);

xij             = xr-coord.x(j);

zij             = zr-coord.z(j);
rij             = sqrt(xij.^2+zij.^2);
drdn            =(xij.*vnx+zij.*vnz)./rij;

t22             = 1i/4*ksi.*besselh(1,2,ksi.*rij).*drdn;
gn22            = G22_SH(ksi,rij,mu);

%integracion gaussiana donde se requiere
ig              = find((rij<=1*para.npplo*drj).*(rij~=0));
t22(ig)         = T22_SH_r_small_ij(coord,indir(ig),j,ksi,gaussian);
gn22(ig)        = G22_SH_r_small_ij(coord,xr(ig),zr(ig),j,ksi,gaussian,mu);

%cuando rij=0
ig0             = find(rij==0);
signo1          = (coord.Xim(indir(ig0),1)==1) - (coord.Xim(indir(ig0),1)==2);
t22(ig0)        = signo1*0.5/drj;
gn22(ig0)       = G22_SH_r_eq_0(ksi,drj,mu);

%fuente imagen
zij             = zr+coord.z(j); 
rij             = sqrt(xij.^2+zij.^2);
drdn            =(xij.*vnx+zij.*vnz)./rij;

t221            = 1i/4*ksi.*besselh(1,2,ksi.*rij).*drdn;
gn221           = G22_SH(ksi,rij,mu);

%integracion gaussiana donde se requiere
ig              = find((rij<=1*para.npplo*drj).*(rij~=0));
coord.z(j)      = -coord.z(j);
coord.vnz(j)    = -coord.vnz(j);
t221(ig)        = T22_SH_r_small_ij(coord,indir(ig),j,ksi,gaussian);
gn221(ig)       = G22_SH_r_small_ij(coord,xr(ig),zr(ig),j,ksi,gaussian,mu);
%cuando rij=0 impossible con la fuente imagen

t22             = (t22 + t221 ).*drj;
gn22            = (gn22+ gn221).*drj;


%  t22(jj)  = ((m==im0) - (m~=im0))*t22(jj);%medio 1 siempre el medio de mas bajo indice