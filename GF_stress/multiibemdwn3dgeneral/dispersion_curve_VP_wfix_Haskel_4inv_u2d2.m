function [vg,vp,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_u2d2(para,wi)

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
err     = 1e-3;


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
    if  pol==1
        kf          = wj/vmin;
    else
        kf          = wj/vmin*1.3;
    end
        
    nk1         = max(round(150/nmodej),30);
%         nk1         = 15;


    nk          = nmodej*nk1;
    nk0         = nk;
    
    while passtest==1
        
        if exist('k20','var')
            %init
            k20(k20<ki) = [];
%             k20(end)    = kf;
%             k20(k20>kf) = [];
            if ~isempty(k20)
                newk2       = zeros(1,length(k20)*nk1);
                
                newk2(1:nk1)= linspace(ki,k20(1),nk1);
                for i=1:length(k20)-1
                    newk2((1:nk1)+i*nk1)=linspace(k20(i),k20(i+1),nk1);
                end
                DWN0.k2 	= unique(newk2);
            else
                DWN0.k2 	= linspace(ki,kf,nk);
            end
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
             if nk/nk0>200
                passtest=0;
                
             end
        elseif  nind>nmodej
            if sum(wj-w0c>=-1e-3)
                %problema de identificacion del primer punto del modo # nmodej
                nmodej=sum(wj-w0c>=-1e-3);
            else
                disp('¿demasiados modos?')
            end

        else
            %mejora la convergencia pero estaria mucha mas sencillo en
            %interpolation lin
            det0  	= real(tmp(indd1));
            det1  	= real(tmp(indd1+1));
            k200  	= DWN0.k2(indd1);
            k201  	= DWN0.k2(indd1+1);
            
            a       = (det1-det0)./(k201-k200).';
            b       =  det1-a.*k201.';
            k20     = -b./a;

            DWN1    = DWN0;
            DWN1.k2	= k20.';
            if  pol==1
                tmp = mode_Love(para,DWN1);
            else
                tmp = mode_Rayleigh_2(para,DWN1);
            end
            tmp(tmp==inf)=0;
            
            %escoger el det de mismo signo y reemplazar
            ind01   = (sign(det0)==sign(real(tmp)));
            det0( ind01)=tmp(ind01);
            det1(~ind01)=tmp(~ind01);
            k200( ind01)=k20(ind01);
            k201(~ind01)=k20(~ind01);
            %reinterpolar
            a       = (det1-det0)./(k201-k200).';
            b       =  det1-a.*k201.';
            k20up   = -b./a;
            
            vp0     = wj./k20;
            vp1     = wj./k20up;
             
            %se recalcula puntos por los cuales el error entre el punto 
            %encontrado y el punto interpolado es mayor al error
            while sum(abs(vp1-vp0)>err)>0
                k20     = k20up;
                DWN1    = DWN0;
                DWN1.k2	= k20.';
                if  pol==1
                    tmp = mode_Love(para,DWN1);
                else
                    tmp = mode_Rayleigh_2(para,DWN1);
                end
                tmp(tmp==inf)=0;
                
                %escoger el det de mismo signo y reemplazar
                ind01   = (sign(det0)==sign(real(tmp)));
                det0( ind01)=tmp(ind01);
                det1(~ind01)=tmp(~ind01);
                k200( ind01)=k20(ind01);
                k201(~ind01)=k20(~ind01);
                %reinterpolar
                a       = (det1-det0)./(k201-k200).';
                b       =  det1-a.*k201.';
                k20up   = -b./a;
                
                vp0         = wj./k20;
                vp1         = wj./k20up;
            end
            k20     = k20up;
            passtest= 0;
        end
        
    end
    k01(j,:)= [fliplr(k20.') zeros(1,nmode-nmodej)].';
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
    if j==1 
        if wi(1)==0
            indi=1;
        end
    end
    kj              = k01(indi:indf,j);
    wj              = wi(indi:indf).';
    plot(kj,wj,'')
    
    ikmax(j)        = length(kj);
    kx(j,1:ikmax(j))= kj;
    vp(j,1:ikmax(j))= wj./kj;
    if j==1 
        if wi(1)==0
            vp(1,1)=para.sub(para.nsubmed).bet;
        end
    end
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
    plot(f1(j,1:ikmax(j)),vp(j,1:ikmax(j)),'r')%
end
xlabel('f');ylabel('vp');

figure(207);hold on;
for j=1:nmode
    plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'c')%
end
xlabel('f');ylabel('vg');