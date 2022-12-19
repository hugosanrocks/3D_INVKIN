function [U,t,DWN] = campo_ref_OP3D_DWN(iinc,coord,para,DWN,polOP,kxk)
%% 3D Incompleto (tomado de PSV)
% Desplazamientos y tracciones para ondas planas en estratificado 3D

% funcion que da los campos de desplazamientos y las tractiones
% para una incidencia de ondas planas en un semi-espacio
% xf,zf coordenadades del origen de la onda (fase nula)
% #r relativo a los receptores
% se incluye de una ves los receptores reales para el calculo de los
% campos incidente

n       = coord.nbpt;
U       = zeros(3,n);
t       = zeros(3,n);
indi  	= coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1

%posicion receptores virtuales (puntos de colocacion) y reales 
nrv     = DWN.nxrv;
xr      = DWN.xr;
yr      = DWN.yr;
zr0     = DWN.zr0;
izr0    = DWN.izr0;
salu    = DWN.salu;
sals    = DWN.sals;

%posicion fuente real
coordf.xs       = para.xs(iinc);
coordf.ys       = para.ys(iinc);
coordf.zs       = para.zs(iinc);

%[u,S]           =  calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,polOP);
kzsigno     = para.kzsigno(iinc);
phi         = para.phi(iinc);
[u,s]     	= calcul_US_DWN_3D_Ncapas_HS_incOP(para,xr,yr,zr0,izr0,salu,sals,coordf,kxk,polOP,kzsigno,phi);
            
% receptores reales
DWN.uxz0(:,:,iinc)   = u(:,nrv+1:end);
% puntos de colocacion
U(:,indi)           = u(:,1:nrv);
s                   = s(:,:,1:nrv);

%calculo de la traccion en los puntos de colocacion solamente
if isfield(coord,'vnx')
vn(1,:)	= coord.vnx(indi);
vn(2,:)	= coord.vny(indi);
vn(3,:)	= coord.vnz(indi);

t(1,indi)  = squeeze(s(1,1,:)).'.*vn(1,:)+squeeze(s(1,2,:)).'.*vn(2,:)+squeeze(s(1,3,:)).'.*vn(3,:);
t(2,indi)  = squeeze(s(1,2,:)).'.*vn(1,:)+squeeze(s(2,2,:)).'.*vn(2,:)+squeeze(s(2,3,:)).'.*vn(3,:);
t(3,indi)  = squeeze(s(1,3,:)).'.*vn(1,:)+squeeze(s(2,3,:)).'.*vn(2,:)+squeeze(s(3,3,:)).'.*vn(3,:);
end
end