function [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv_adapt(para,wi)

% funcion que busqua a establecer las curvas de dispersion de velocidad de grupo
% por un medio estratificado sobre un semi espacio
% se busca por cada w(wfix) 
% los ceros de la fase de un sub-determinante de la matriz de haskel 
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
%
% diferencia con dispersion_curve_wfix_Haskel_4inv:
% se da un criterio de convergencia para estar seguro que la velocidad de
% grupo esta bien calculada

pol     = para.pol;

beta    = zeros(para.nsubmed,1);
alpha   = zeros(para.nsubmed,1);
if isfield(para,'reg')
    para.sub= para.reg(1).sub;
end
for ms=1:para.nsubmed
    beta(ms)  = para.sub(ms).bet;
    alpha(ms) = para.sub(ms).alpha;
    para.sub(ms).C(6,6)  = para.sub(ms).rho*para.sub(ms).bet^2;
end

vmin    = min(beta);
vmax    = max(beta);

vmax    = (1-1e-8)*vmax;

wi      = unique(wi);
dwi     = min(diff(wi));
dwi     = max(dwi*1e-3,1e-2);
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

%factor que permite encontrar mas o menos rapidamente las raices con la precision deseada
%depende de la maquina y es por la programacion vectorial de matlab
%para pasar a fortran, ocupar nk=3;
nk      = 150;
nk0     = nk;

indk    = 1:(nk-1);
indk0   = indk;

j=1;
while j<=nf
    wj          = wig(j);
    nmodej      = sum(w0c<=wj);
    if j==1
        k2          = [wj/vmax, min(wj/vmin,wmax)*1.2];
    else
        k200(k200<(wj/vmax))=[];
        k200        = unique(k200);
        k2          = [wj/vmax, k200,min(wj/vmin,wmax)*1.2];
    end
    DWN0        = struct('omegac',wj);
    ik          = 1;
    k20         = zeros(1,nmodej);
    i           = 1;
    while i<=(length(k2)-1)
        DK      = (k2(i+1)-k2(i));
        k2i     = k2(i);
        k2f     = k2(i+1);
        a       = 1;
        force   = 0;
        
        while ((DK>1e-3*dwi) || (a==1)) && force==0
            a       = 0;%pour rentrer ds la boucle au moins une fois
            k2p     = linspace(k2i,k2f,nk);
            DWN0.k2	= k2p;
%             det0=0;
            if  pol==1
                det0=mode_Love(para,DWN0);
            else
                det0=mode_Rayleigh_2(para,DWN0);
            end
            det0(det0==inf)=0;
            det_ph    = angle(det0)-pi/2;
            indd1   = logical(det_ph(1:(nk-1)).*det_ph(2:nk)<=0);
            indd1   = indk(indd1);
            if isempty(indd1)
                break
            elseif length(indd1)==2
                indd1=indd1(1);
                k2i     = k2p(indd1);
                k2f     = k2p(indd1+1);
                k2(i)   = k2f;
                i       = i-1;
            elseif length(indd1)>2
                %disp('pb of precision')
                force   = 1;
                indd1   = indd1(1);
            else
                k2i     = k2p(indd1);
                k2f     = k2p(indd1+1);
            end
            DK      = (k2f-k2i);
        end
       if ~isempty(indd1)
            %interp quad
            k20(ik)	= cherche_zero(DWN0.k2,real(det0).',indd1).';
            ik      = ik+1;
        end
        i       = i+1;
    end
    
    nind    = nmodej-sum(k20==0)-sum(isnan(k20));
    if wj==0
        k20	= 0;
        nind= 1;
    end
    if nind<nmodej 
        nk      = 5*nk;
        indk    = 1:(nk-1);
        if nk>100*nk0
            break
        end
    elseif length(k20)>nmodej
        nk      = nk0;
        indk    = indk0;
%         k01(j,:)= [flipud(k20(2:end).');zeros(nmode-nmodej,1)].';
%         j       = j+1;
%         k200    = k20(2:end);
    else
        nk      = nk0;
        indk    = indk0;
        k01(j,:)= [flipud(k20.');zeros(nmode-nind,1)].';
        j       = j+1;
        k200    = k20;
    end
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
%     plot(kj,wj,'k')
% 
%     figure(207);hold on;
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'r')
end