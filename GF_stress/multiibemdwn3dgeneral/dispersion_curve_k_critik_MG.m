function [k20,w0]=dispersion_curve_k_critik_MG(para,wmax,v)

% funcion que busqua los arrenques sobre la pendiente w=k*vmax
% se busca los zeros del determinante de la matriz del DWN por cada k

fj          = linspace(0,wmax/2/pi,1e3);
k20         = 2*pi*fj/v;
DWN0        = struct('k2',k20);
DWN0.omegac	= 2*pi*fj;

for ic=1:para.nsubmed
    para.reg(1).sub(ic).ksi  = DWN0.omegac/para.reg(1).sub(ic).bet;
end

if  para.pol==1
    DWN0     = calcul_A_DWN_SH_Ncapas_HS(para,DWN0);
else
    for ic=1:para.nsubmed
        para.reg(1).sub(ic).kpi  = DWN0.omegac/para.reg(1).sub(ic).alpha;
    end
    DWN0     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN0);
end

nw=length(DWN0.omegac);
tmp=zeros(nw,1);
for i=1:nw
    tmp(i)=(det(DWN0.A_DWN(:,:,i)));
end

%         figure(203);hold on;plot(DWN0.omegac,real(tmp),'k');plot(DWN0.omegac,imag(tmp),'r');plot(DWN0.omegac,abs(tmp),'')

w0  = cherchermin_w_k2V(abs(tmp),para,DWN0,v,100);
w0  = [0;w0];
k20	= w0/v;
% figure(205);hold on;plot(w0/vmax,w0,'.r','markersize',6)