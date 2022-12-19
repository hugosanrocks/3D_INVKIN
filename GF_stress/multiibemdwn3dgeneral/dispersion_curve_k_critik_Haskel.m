function [k20,w0]=dispersion_curve_k_critik_Haskel(para,wmax,v,nw)

% funcion que busqua los arrenques sobre la pendiente w=k*vmax
% se busca los ceros de la fase de un sub-determinante de la matriz de haskel
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
% los arranques son bastante regulares

if ~exist('nw','var'); nw = 1e3; end
wj          = linspace(0,wmax,nw);
k20         = wj/v;

DWN0.k2 	= k20;
DWN0.omegac	= wj;

if  para.pol==1
    det0= mode_Love(para,DWN0);
else
    det0= mode_Rayleigh_2(para,DWN0);
end
det_ph	= angle(det0)-pi/2;
indd1   = logical(det_ph(1:(nw-1)).*det_ph(2:nw)<=0);
indk    = 1:(nw-1);
indd1   = indk(indd1);
if ~isempty(indd1)
    [k20,~]     = cherche_zero(DWN0.k2,real(det0).',indd1);
    k20         = k20.';
    if indd1(1)==1
        k20(1)  = 0;
    else
        k20=[0;k20];
    end
else
    k20  = 0;
end
w0      = k20*v;