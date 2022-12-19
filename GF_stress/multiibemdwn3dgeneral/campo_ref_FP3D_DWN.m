function [U,t,DWN]=campo_ref_FP3D_DWN(iinc,coord,para,DWN)
% funcion que da los campos de desplazamientos y las tractiones
% para una fuente vectorial en un espacio estratificado
% se incluye de una ves los receptores reales para el calculo de los
% campos incidente

n       = coord.nbpt;
t       = zeros(3,n);
U       = zeros(3,n);

indi  	= coord.indm(1,1).ind;%true-false de los puntos perteneciendo a m=1

%posicion receptores virtuales (puntos de colocacion) y reales 
nrv             = DWN.nxrv;

salu            = DWN.salu;
sals            = DWN.sals;

rec.xr          = DWN.xr;
rec.yr          = DWN.yr;
rec.zr          = DWN.zr;
rec.zr0         = DWN.zr0;
rec.izr0        = DWN.izr0;
rec.zricrall    = DWN.zricrall;
rec.icrall      = DWN.icrall;
nrec            = length(DWN.xr);

%fuente real
coordf.xs       = para.xs(iinc);
coordf.ys       = para.ys(iinc);
coordf.zs       = para.zs(iinc);
fij             = para.fij(iinc,:);

[u_f1,s_f1,u_f2,s_f2,u_f3,s_f3] = calcul_US_DWN_3D_polar_Ncapas_HS2(para,rec,salu,sals,coordf,fij,DWN);

u = u_f1+u_f2+u_f3;
S = s_f1+s_f2+s_f3;

% receptores reales
DWN.U0(:,:,iinc)	= u(:,nrv+1:end);
if length(sals)>nrv && sals(nrv+1)
    S1              = squeeze(reshape(S(:,:,nrv+1:nrec),1,9,nrec-nrv));
    S1              = S1(:,DWN.inds);
    DWN.S0(:,:,iinc)= S1;
end

% puntos de colocacion
U(:,indi)           = u(:,1:nrv);
S                   = S(:,:,1:nrv);

%calculo de la traccion en los puntos de colocacion solamente
vn(1,:)	= coord.vnx(indi);
vn(2,:)	= coord.vny(indi);
vn(3,:)	= coord.vnz(indi);

t(1,indi)  = squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:)+squeeze(S(1,3,:)).'.*vn(3,:);
t(2,indi)  = squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:)+squeeze(S(2,3,:)).'.*vn(3,:);
t(3,indi)  = squeeze(S(1,3,:)).'.*vn(1,:)+squeeze(S(3,2,:)).'.*vn(2,:)+squeeze(S(3,3,:)).'.*vn(3,:);