function [vg,f1,ikmax]=dispersion_curve_wfix_Haskel_4inv_adapt2(para,wi) %#codegen

% funcion que busqua a establecer las curvas de dispersion de velocidad de grupo
% por un medio estratificado sobre un semi espacio
% se busca por cada w(wfix)
% los ceros de la fase de un sub-determinante de la matriz de haskel
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
%
% diferencia con dispersion_curve_wfix_Haskel_4inv:
% se trata de una programacion no orientada matlab, cuando es compilado
% tendria que ser mas eficiente ya que se busca a converger mas rapidamente
% a la solucion con bissection ocupando las soluciones al paso anterior
% evita asi calcular un gran numero de puntos no deseado
% se cambia las entradas para que no haya estructuras

% reescrito para mex

coder.varsize('k200','kj','tmp','indk','k2p');


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
dwi     = max(dwi*1e-3,1e-4);
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
%para pasar a fortran, ocupar nk=2;
nk      = 3;
nk0     = nk;

indk    = 1:(nk-1);
indk0   = indk;
k200    = zeros(1,nmode);
j=1;
while j<=nf
    wj          = wig(j);
    nmodej      = sum(w0c<=wj);
    if j==1
        k2          = [wj/vmax, min(wj/vmin,wmax)*1.2];
    else
        k200    = unique(k200);
        k200(k200<(wj/vmax))=[];
        
        k2      = [wj/vmax, k200,min(wj/vmin,wmax)*1.2];
    end
    DWN0        = struct('omegac',wj,'k2',0);
    ik          = 1;
    k20         = zeros(1,nmodej);
    i           = 1;
    while i<=(length(k2)-1)
        DK      = (k2(i+1)-k2(i));
        k2i     = k2(i);
        k2f     = k2(i+1);
        a       = 1;
        force   = 0;
        indd1   = 0;
        det0    = 0;
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
                k2i     = k2p(indd1(1));
                k2f     = k2p(indd1(1)+1);
                k2(i)   = k2f;
                i       = i-1;
            elseif length(indd1)>2
                %disp('pb of precision')
                force   = 1;
                indd1   = indd1(1);
            else
                k2i     = k2p(indd1(1));
                k2f     = k2p(indd1(1)+1);
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
    kj              = k01(indi(1):indf(1),j);
    if mod(length(kj),2)~=0
        indi(1)        = indi(1)+1;
        kj         	= k01(indi(1):indf(1),j);
    end
    wj0              = wig(indi(1):indf(1)).';
    ikmax(j)        = length(kj)/2;
    tmp             = diff(wj0)./diff(kj);
    vg(j,1:ikmax(j))= tmp(1:2:end);
    
    f1(j,1:ikmax(j))= (wj0(1:2:(end-1))+wj0(2:2:end))/2/2/pi;
end