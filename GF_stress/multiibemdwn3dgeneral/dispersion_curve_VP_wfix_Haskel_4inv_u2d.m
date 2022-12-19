function [vg,vp,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_u2d(para,wi)

% funcion que busqua a establecer las curvas de dispersion de velocidad de
% fase
% por un medio estratificado sobre un semi espacio
% se busca por cada w(wfix)
% los ceros de la fase de un sub-determinante de la matriz de haskel
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
%
% diferencia con dispersion_curve_wfix_Haskel_4inv:
% no mas las velocidad de fase
err     = 1e-2;


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
wmax    = max(wi);
nf      = length(wi);
%&precalculo lista indices
if para.pol==2
    jj=ind4subdetab(para.nsubmed);
    para.jj=jj;
end
[k2c,w0c]   = dispersion_curve_k_critik_Haskel(para,wmax,vmax);
nmode       = length(k2c);
k01         = zeros(nf,nmode);

for j=nf:-1:1
    wj          = wi(j);
    if wj==0
        k20     = 0;
        %automaticamente hecho porque k20=0
        % nmodej  = 1;
        % k01(j,:)= [flipud(k20);zeros(nmode-nmodej,1)].';
        continue;
    end
    nmodej      = sum(w0c<=wj);
    DWN0        = struct('omegac',wj);
    
    passtest    = 1;
    %limites en k dados por asintotas
    ki          = wj/vmax;
    kf          = wj/vmin*1.3;
    
    nk1         = max(round(150/nmodej),15);
%         nk1         = 15;


    nk          = nmodej*nk1;
    nk0         = nk;
    
    while passtest==1
        
        if exist('k20','var')
            %init
            k20(k20<ki)=[];
            k20(k20>kf)=[];

            newk2       = zeros(1,length(k20)*nk1);
            
            newk2(1:nk1)= linspace(ki,k20(1),nk1);
            for i=1:length(k20)-1
                newk2((1:nk1)+i*nk1)=linspace(k20(i),k20(i+1),nk1);
            end
            DWN0.k2 	= unique(newk2);
        else
            DWN0.k2 	= linspace(ki,kf,nk);
        end
        nk      = length(DWN0.k2);
        indk    = 1:(nk-1);
        
        if  pol==1
            tmp = mode_Love(para,DWN0);
        else
            tmp = mode_Rayleigh_2(para,DWN0);
        end
        tmp(tmp==inf)=0;
        tmp1    = angle(tmp)-pi/2;
        indd1   = logical(tmp1(1:(nk-1)).*tmp1(2:nk)<=0);
        indd1   = indk(indd1);
        nind    = length(indd1);
        if nind<nmodej
            nk1     = nk1*10;
            nk    	= nmodej*nk1;
        else
            %mejora la convergencia pero estaria mucha mas sencillo en
            %programacion lineal
            [k20,k20l] 	= cherche_zero(DWN0.k2,real(tmp).',indd1);
            vp0         = wj./DWN0.k2(indd1);            
            vp1         = wj./k20;
            
            %se recalcula puntos por los cuales el error entre el punto mas
            %encontrado y el punto interpolado es mayor al error
            ind0        = find(abs(vp1-vp0)>err);
            %se recalcula puntos por los cuales hay poco puntos entre
            %varios ceros -> mala interpolacion
            pbex        = find(diff(indd1)<4);
            pbex        = unique([pbex, pbex+1]);
            ind0        = unique([pbex, ind0]);
            while ~isempty(ind0)
                nk1     = max(round(150/length(ind0)),10);
                nk2     = nk1*length(ind0);
                indk1   = 1:(nk2-1);
                
                newk2   = zeros(1,nk2);
                for i=1:length(ind0)
                    newk2((1:nk1)+(i-1)*nk1)=linspace(DWN0.k2(indd1(ind0(i))),DWN0.k2(indd1(ind0(i))+1),nk1);
                end
                %cuidado acabamos de introducir intervalos grandes entre
                %cada bloque de nk1. En esos intervalos se puede encontrar
                %zeros que ya tenemos bien. Hay que excluir estos indices
                %(vea despues) la lista de los indices son:
                index   = nk1:nk1:nk2;
                
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
                for iex=1:length(index)
                    indd1(indd1==index(iex))=[];
                end

%                 vp0(ind0)	= wj./newk2(indd1);
                k21         = zeros(1,length(ind0));
                k21l        = zeros(1,length(ind0));
                for i0=1:length(ind0)
                    [k21(1,i0),k21l(1,i0)]   	= cherche_zero(DWN0.k2((1:nk1)+(i0-1)*nk1),real(tmp((1:nk1)+(i0-1)*nk1)).',indd1(i0)-(i0-1)*nk1);
                end
                
                k20(ind0)   = k21;
                k20l(ind0)  = k21l;

                vp0         = wj./k20l;
                vp1         = wj./k20;
                ind00       = find(abs(vp1(ind0)-vp0(ind0))>err);
                ind0        = ind0(ind00);%no se debe cambiar el indice que es el del vector inicial
                indd1(ind0) = indd1(ind00);%se cambia para el uso en el ciclo de arriba
            end
            passtest= 0;
        end
        
    end
    k01(j,:)= [fliplr(k20) zeros(1,nmode-nmodej)].';
end

vp      = zeros(nmode,nf);
vg      = zeros(nmode,nf);
f1      = zeros(nmode,nf);
kx      = zeros(nmode,nf);
ikmax   = zeros(nmode,1);

figure(205);hold on;

for j=1:nmode
    indi=find(k01(:,j)~=0,1,'first');
    indf=find(k01(:,j)~=0,1,'last');
    
    kj              = k01(indi:indf,j);
    wj              = wi(indi:indf).';
    plot(kj,wj,'')
    
    ikmax(j)        = length(kj);
    kx(j,1:ikmax(j))= kj;
    vp(j,1:ikmax(j))= wj./kj;
    f1(j,1:ikmax(j))= wj/2/pi;
    
    vg(j,1:ikmax(j))=  gradient(wj,kj);
end

wj=linspace(0,max(wi),3);
for j=1:para.nsubmed
    plot(wj/para.sub(j).bet,wj,'r')
end

xlabel('k');ylabel('w');

figure(206);hold on;
for j=1:nmode
    plot(vp(j,1:ikmax(j)),'r')%f1(j,1:ikmax(j)),
end
xlabel('f');ylabel('vp');

figure(207);hold on;
for j=1:nmode
    plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'r')%
end
xlabel('f');ylabel('vg');