function [U,t,DWN] = campo_ref_PSV_OP_DWN(iinc,coord,para,DWN,polOP,kxk)

% funcion que da los campos de desplazamientos y las tractiones
% para una incidencia de ondas planas PSV en un semi-espacio
% xf,zf coordenadades del origen de la onda (fase nula)
% #r relativo a los receptores
% se incluye de una ves los receptores reales para el calculo de los
% campos incidente

n       = coord.nbpt;
t       = zeros(2,n);
U       = zeros(2,n);

indi  	= coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1

%posicion receptores virtuales (puntos de colocacion) y reales 
nrv     = DWN.nxrv;
xr      = DWN.xr;
zr0     = DWN.zr0;
izr0    = DWN.izr0;
salu    = DWN.salu;
sals    = DWN.sals;

%posicion fuente real
coordf.xs       = para.xs(iinc);
coordf.zs       = para.zs(iinc);

[u,S]           =  calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,polOP);
% receptores reales
DWN.uxz0(:,:,iinc)   = u(:,nrv+1:end);
% puntos de colocacion
U(:,indi)           = u(:,1:nrv);
S                   = S(:,:,1:nrv);

%calculo de la traccion en los puntos de colocacion solamente
vn(1,:)	= coord.vnx(indi);
vn(2,:)	= coord.vnz(indi);

t(1,indi)= squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:);
t(2,indi)= squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:);
