% function dispersion_curve_wfix_MG(para)

% funcion que busqua a establecer las curvas de dispersion por un medio
% estratificado sobre un semi espacio
% se busca los zeros del determinante de la matriz del DWN por cada w(wfix)
% en un gran intervalo busquando un cambio de signo del determinante
% (los acercamientos a zero sin cambiar de signo no estan considerados)
% comentarios: poco eficiente, puede pasar a lado de puntos
% pb debido a que los modos de volumen se encuentren tambien y hay cruce: 
% 2 ceros muy cercanos que no se ven

tic
wmax    = 2*pi*para.fmax;
para    = Vp2Cij(para);
pol     = para.pol;
ncapas  = para.nsubmed;

beta    = zeros(ncapas,1);
for ms=1:ncapas
    beta(ms)  = para.reg(1).sub(ms).bet;
    para.reg(1).sub(ms).Ci 	= para.reg(1).sub(ms).C;
end
vmin    = min(beta);
vmax    = max(beta)*(1+1e-6);
vdeep   = beta(end)*(1-1e-6);

%numero de frecuencias a buscar los k2 soluciones
nf      = 200; 
dw      = 2*pi*para.fmax/nf;

[k2c,wc]= dispersion_curve_k_critik_MG(para,wmax,vdeep);
nmode 	= length(k2c)+pol*ncapas;
k2sol	= zeros(nf,nmode);
nk      = 1e2*nmode;


for j=1:nf
    wj          = j*dw + 0.0001*dw*(j==0);
    DWN         = struct('omegac',wj);
    DWN.k2    	= linspace(wj/vmax,min(wj/vmin,wmax)*1.2,nk);
    
    for ic=1:ncapas
        para.reg(1).sub(ic).ksi  = DWN.omegac/para.reg(1).sub(ic).bet;
    end
    if  pol==1
        DWN = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
    else
        for ic=1:ncapas
            para.reg(1).sub(ic).kpi  = DWN.omegac/para.reg(1).sub(ic).alpha;
        end
        DWN = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
    end
    
    tmp = zeros(nk,1);
    for i=1:nk
        tmp(i)=(det(DWN.A_DWN(:,:,i)));
    end
    
    k20         = cherchermin_w(abs(tmp),para,DWN,100);
    k20(k20==0) = [];
    nind        = length(k20);
    k2sol(j,:)    = [flipud(k20);zeros(nmode-nind,1)];%k20 crepe pour continuite des modes
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(205);hold on;
plot(k2c,wc,'.r','markersize',6)
j       = 1:nf;
wj      = j*dw + 0.0001*dw*(j==0);
for j=1:nf
    iind=find(k2sol(j,:)~=0,1,'last');
    plot(k2sol(j,1:iind),wj(j)*ones(1,iind),'.r','markersize',6)
end
xlabel('k_x');
ylabel('\omega');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica asintotas                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2  = 0:1e3;
for ic=1:ncapas
    C   = para.reg(1).sub(ic).C;
    rho = para.reg(1).sub(ic).rho;
    vT  = sqrt(C(6,6)/rho);
    vL  = sqrt(C(1,1)/rho);
    if  para.pol==1
        w=k2*vT;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2(1:i0),w(1:i0),'k');
    else
        w=k2*vT;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2(1:i0),w(1:i0),'');
        w   = k2*vL;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2(1:i0),w(1:i0),'r');
    end
end
