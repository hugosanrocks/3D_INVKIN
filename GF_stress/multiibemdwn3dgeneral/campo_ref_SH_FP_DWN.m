function [U,t,DWN]=campo_ref_SH_FP_DWN(iinc,coord,para,DWN)
% funcion que da los campos de desplazamientos y las tractiones
% para una fuente puntual SH en un multi-estarto
% xf,zf coordenadades del origen de la onda (fase nula)
% #r relativo a los receptores
% se incluye de una ves los receptores reales para el calculo de los
% campos incidente

n       = coord.nbpt;
t       = zeros(1,n);
U       = zeros(1,n);

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

[u,S]           = calcul_US_DWN_SH_Ncapas_HS(para,xr,zr0,izr0,salu,sals,coordf,DWN);
DWN.uy0(:,iinc) = u(nrv+1:end);
DWN.s0(:,:,iinc)= S(DWN.inds,nrv+1:end);
U(indi)         = u(1:nrv);
S               = S(:,1:nrv);
% S = [SxyKW,SzyKW](xr,zr,w)

%calculo de la traccion en los puntos de colocacion solamente
vn(1,:)	= coord.vnx(indi);
vn(2,:)	= coord.vnz(indi);

t(indi)  = S(1,:).*vn(1,:)+S(2,:).*vn(2,:);
