function dispersion_curve_kfix_MG(para)

% funcion que busqua a establecer las curvas de dispersion por un medio
% estratificado sobre un semi espacio
% se busca los zeros del determinante de la matriz del DWN (matriz global)
% por cada k(kfix) en un gran intervalo 
% busquando un cambio de signo del determinante
% (los acercamientos a zero sin cambiar de signo no estan considerados)
% comentarios: poco eficiente, puede pasar a lado de puntos
%              pero nunca ve errores, estable
% la paralelizacion es posible pero no interesante

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
vmax    = max(beta)*(1-1e-6);
vdeep   = beta(end)*(1-1e-6);

%numero de k2 a buscar los w soluciones
nk          = 200;
kmax        = wmax/min(beta);
DK          = kmax/nk;


[k2c,wc]   = dispersion_curve_k_critik_MG(para,wmax,vdeep);
nmode       = length(k2c)+pol*ncapas;
wsol        = zeros(nk,nmode);
nf          = 1e2*nmode;

for j=1:nk
    k20         = j*DK;
    DWN         = struct('k2',k20);
    
    DWN.omegac	= 2*pi*linspace(k20*vmin*.5,min(k20*vmax,wmax),nf+1)/2/pi;
    
    for ic=1:ncapas
        para.reg(1).sub(ic).ksi  = DWN.omegac/para.reg(1).sub(ic).bet;
    end
    if  pol==1
        DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
    else
        for ic=1:ncapas
            para.reg(1).sub(ic).kpi  = DWN.omegac/para.reg(1).sub(ic).alpha;
        end
        DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
    end
    % nw=length(DWN.omegac);tmp=zeros(nw,1);for i=1:nw;tmp(i)=(det(DWN.A_DWN(:,:,i))); end;
    % figure(206);hold on;plot(DWN.omegac,real(tmp),'');plot(DWN.omegac,imag(tmp),'r')

    tmp=zeros(nf+1,1);
    for i=1:nf+1
        tmp(i)=(det(DWN.A_DWN(:,:,i)));
    end
    
    % tmp=determinant_vec_de_mat(DWN.A_DWN);
    % figure(203);hold on;plot(DWN.omegac,real(tmp)/max(abs(real(tmp))),'k');plot(DWN.omegac,imag(tmp)/max(abs(imag(tmp))),'r');plot(DWN.omegac,abs(tmp)/max(abs(tmp)),'');plot(DWN.omegac,angle(tmp)/max(abs(angle(tmp))),'c')
    
    w0          = cherchermin_k(abs(tmp),para,DWN,100);
    w0(w0==0)   = [];
    nind        = length(w0);
    wsol(j,:)    = [w0;zeros(nmode-nind,1)];
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2  = (1:nk)*DK;
figure(205);hold on;
plot(k2c,wc,'.k','markersize',6)

for j=1:nk
    iind=find(wsol(j,:)~=0,1,'last');
    plot(k2(j)*ones(iind,1),wsol(j,1:iind),'.k')
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


% for ic=1:ncapas
%     if  para.pol==1
%         if ic<ncapas
%             %             f0  = vT/para.reg(1).sub(ic).h/4;
%             %             f   =(0:11)*f0;
%             %             figure(205);hold on;plot(0,f*2*pi,'.c');
%             for n=0:10
%                 f   = vT*sqrt((1/para.reg(1).sub(ic).h/4*(n))^2+(k2/2/pi).^2);
%                 figure(205);hold on;plot(k2,f*2*pi,'c');
%             end
%         else
%         end
%     else
%         if ic==1<ncapas
%             f0  = vT/para.reg(1).sub(ic).h/4;
%             f   =(0:11)*f0;
%             figure(205);hold on;plot(0,f*2*pi,'.c');
%             for n=0:10
%                 f   = vT*sqrt((1/para.reg(1).sub(ic).h/4*(2*n+1))^2+(k2/2/pi).^2);
%                 figure(205);hold on;plot(k2,f*2*pi,'c');
%             end
%         end
%     end
% end


%1 capa/Semi-espacio
% for imode =1:5
%     k2  = vT1*imode*pi/(h1*sqrt(vT2^2-vT1^2));
%     w   = k2*vT2;
%     plot(k2,w,'.r')
% end



% %2 capas/Semi-espacio
% h1      = para.reg(1).sub(1).h;
% C       = para.reg(1).sub(1).C;
% rho     = para.reg(1).sub(1).rho;
% vT1     = sqrt(C(6,6)/rho);
% h2      = para.reg(1).sub(2).h;
% C       = para.reg(1).sub(2).C;
% rho     = para.reg(1).sub(2).rho;
% vT2     = sqrt(C(6,6)/rho);
% C       = para.reg(1).sub(3).C;
% rho     = para.reg(1).sub(3).rho;
% vT3     = sqrt(C(6,6)/rho);
%
%
% alpha=(vT2^2/(vT3^2-vT2^2)/vT1^2*(vT3^2-vT1^2));
% for imode2=0:3
%     for imode1=1:3
%         k2  = vT2/(sqrt(vT3^2-vT2^2))*pi/(h1+h2)*(imode1*sqrt(alpha)/6+imode2/4);
%         w   = k2*vT3;
%         plot(k2,w,'oc')
%     end
% end