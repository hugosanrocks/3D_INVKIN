function [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv(para,wi)

% funcion que busqua a establecer las curvas de dispersion de velocidad de grupo
% por un medio estratificado sobre un semi espacio
% se busca por cada w(wfix)
% los ceros de la fase de un sub-determinante de la matriz de haskel
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
%
% diferencia con dispersion_curve_wfix_Haskel:
% indicada para calcular velocidad de grupo en las frecuencias deseada
% se calcula el numero de modos que se espera para cada frecuencia
% si no se encuentra el numero de modo esperado se aumenta el numero de
% punto para la busqueda de cambio de signo
% se considera un numero de punto de tal manera que la velocidad de grupo
% sea calculada con una precision deseada

para.sub= para.reg(1).sub;
pol     = para.pol;

beta    = zeros(para.nsubmed,1);
alpha   = zeros(para.nsubmed,1);
for ms=1:para.nsubmed
    para.sub(ms).C(6,6)  = para.sub(ms).rho*para.sub(ms).bet^2;
    beta(ms)  = para.sub(ms).bet;
    alpha(ms) = para.sub(ms).alpha;
end

vmin    = min(beta);
vmax    = max(beta);

vmax    = (1-1e-8)*vmax;

wi      = unique(wi);
dwi     = min(diff(wi))/10;
% dwi     = max(dwi*1e-3,1e-4);
wig     = zeros(1,2*length(wi));
wig(1:2:end)= wi-dwi;
wig(2:2:end)= wi+dwi;
if max(wig<0)==1
    wig(1)=0;
    wig(2)=1e-3;
end
wmax    = max(wig);

nf      = length(wig);
[k2c,w0c]   = dispersion_curve_k_critik_Haskel(para,wmax,vmax);
nmode       = length(k2c);
k01         = zeros(nf,nmode);
for j=1:nf
% parfor j=1:nf
    wj          = wig(j);
    if wj==0
        k20=0;
        continue;
    end
    nmodej      = sum(w0c<=wj);
    DWN0        = struct('omegac',wj);
    
    passtest    = 1;
    ki          = wj/vmax;
    kf          = min(wj/vmin,wmax)*1.2;
    nk          = nmodej*10;%min(floor(max(dk0/dwi*1e3,nmodej*3)),5e3);
    nk0         = nk;
    while passtest==1
        indk        = 1:(nk-1);
        DWN0.k2 	= linspace(ki,kf,nk);
        if  pol==1
            tmp=mode_Love(para,DWN0);
        else
            tmp=mode_Rayleigh_2(para,DWN0);
        end
        tmp(tmp==inf)=0;
        tmp1    = angle(tmp)-pi/2;
        indd1   = logical(tmp1(1:(nk-1)).*tmp1(2:nk)<=0);
        indd1   = indk(indd1);
        nind    = length(indd1);
        if nind<nmodej
            nk      = nk*10;
            if nk/nk0>200
                %                 k20=[k2c(nmodej);k20];
                passtest=0;
            end
        elseif nind>nmodej
            %k20=k20(1:nmodej);
            passtest=0;
            disp('pb')
        else
            nk1     = max(round(150/nmodej),10);
            nk2     = nk1*nmodej;
            dk      = DWN0.k2(2)-DWN0.k2(1);
            indk1        = 1:(nk2-1);
            
            while dk>dwi*1e-3/vmin
                newk2=zeros(1,nk2);
                for i=1:nmodej
                    newk2((1:nk1)+(i-1)*nk1)=linspace(DWN0.k2(indd1(i)),DWN0.k2(indd1(i)+1),nk1);
                end
                DWN0.k2=newk2;
                if  pol==1
                    tmp=mode_Love(para,DWN0);
                else
                    tmp=mode_Rayleigh_2(para,DWN0);
                end
                tmp(tmp==inf)=0;
                tmp1    = angle(tmp)-pi/2;
                indd1   = logical(tmp1(1:(nk2-1)).*tmp1(2:nk2)<=0);
                indd1   = indk1(indd1);
                dk      = DWN0.k2(2)-DWN0.k2(1);
            end
%             k20   	= DWN0.k2(indd1).';
            k20   	= cherche_zero(DWN0.k2,real(tmp).',indd1).';
            passtest=0;
        end
        
    end
    
    k01(j,:)= [flipud(k20);zeros(nmode-nmodej,1)].';
end
vg      = zeros(nmode,nf);
f1      = zeros(nmode,nf);
ikmax   = zeros(nmode,1);

for j=1:nmode
    indi=find(k01(:,j)~=0,1,'first');
    indf=find(k01(:,j)~=0,1,'last');
    
    kj              = k01(indi:indf,j);
    if mod(length(kj),2)~=0
        indi        = indi+1;
        kj         	= k01(indi:indf,j);
    end
    wj              = wig(indi:indf).';
    ikmax(j)        = length(kj)/2;
    tmp             = diff(wj)./diff(kj);
    vg(j,1:ikmax(j))= tmp(1:2:end);
    
    f1(j,1:ikmax(j))= (wj(1:2:(end-1))+wj(2:2:end))/2/2/pi;
%     figure(205);hold on;
%     plot(kj,wj,'')
%     figure(206);hold on;
%     plot(wj/2/pi,wj./kj,'')
%     figure(207);hold on;
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'r')
end
% toto=0;