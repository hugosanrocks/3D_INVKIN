function [vg,vp,f1,ikmax,kx]=dispersion_curve_VP_wfix_Haskel_4inv_d2u3(para,wi)

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

vmin    = (1+1e-8)*min(beta);
vmax    = (1-1e-8)*max(beta);

wi      = unique(wi);
wmax    = max(wi);
nf      = length(wi);

test    = 0;
nw      = 1e4;
while test==0
    [k2c,w0c]   = dispersion_curve_k_critik_Haskel(para,wmax,vmax,nw);
    
    nmode       = length(k2c);
    if nw>nmode*20
        test    = 1;
    else
        nw      = nmode*20*2;
    end
end
k01         = zeros(nf,nmode);
pbw         = zeros(nf,1);

for j=1:nf
%     disp(j)
    wj          = wi(j);
    if j<nf
        dw          = wi(j)-wi(j+1);
    end
    if wj==0
        k20     = 0;
        %automaticamente hecho porque k20=0
        % nmodej  = 1;
        % k01(j,:)= [flipud(k20);zeros(nmode-nmodej,1)].';
        continue;
    end
    nmodej      = sum(w0c<=wj);
    DWN0        = struct('omegac',wj);
    
    %limites en k dados por asintotas
    ki          = wj/vmax;
    if  pol==1
        kf    	= wj/vmin;
    else
        kf  	= wj/vmin*1.4;
    end
    
    nk1         = max(round(150/nmodej),15);
    nk          = nmodej*nk1;
    nk0         = nk;
    
%     %primer intento
%     %extrapolacion de los valores anteriores para prever mejor los limites
%     %de busqueda
%     if j>3
%         nmodejp1= find(k01(j-1,:)~=0,1,'last');
%         k0g     = interp1(wi(j-3:j-1),k01(j-3:j-1,1:nmodejp1),wi(j),'spline','extrap');
%         dk0g     = diff(fliplr(k0g));
%         %cuando los intervalos son chiquitos se amplia el ancho para no
%         %perder los modos
%         for i=1:length(dk0g)-1
%             if dk0g(i)<dk0g(i+1)/(nk1*2)
%                 %el intervalo corriente es muy pequeno
%                 %se amplia
%                 if i>1
%                     if k0g(i)-dk0g(i)>k0g(i-1)
%                         k0g(i)=k0g(i)-dk0g(i);
%                     end
%                 end
%                 k0g(i+1)=k0g(i+1)+dk0g(i);
%             end
%         end
%         if length(k20)~=length(k0g)
%             toto=0;
%         end
%         k20=k0g;
%     end
    if exist('k20','var')
        %init
        k20(k20<ki) = [];

        if ~isempty(k20)
            newk2       = zeros(1,(length(k20)+1)*nk1);
            
            newk2(1:nk1)= linspace(ki,k20(1),nk1);
            for i=1:length(k20)-1
                newk2((1:nk1)+i*nk1)=linspace(k20(i),k20(i+1),nk1);
            end
            i=length(k20);
            newk2((1:nk1)+i*nk1)=linspace(k20(i),kf,nk1);

            DWN0.k2 	= unique(newk2);
        else
            DWN0.k2 	= linspace(ki,kf,nk);
            %se puede mejorar tomando en cuenta el numero de modo entre las
            %asintotas y fijar otras densidades de puntos
        end
    else
        DWN0.k2 	= linspace(ki,kf,nk);
    end
    nk      = length(DWN0.k2);
    indk    = 1:(nk-1);
    
    if  pol==1
        det = mode_Love(para,DWN0);
    else
        det = mode_Rayleigh_2(para,DWN0);
    end
    det(det==inf)=0;
    detf    = angle(det)-pi/2;
    indd1   = logical(detf(1:(nk-1)).*detf(2:nk)<=0);
    indd1   = indk(indd1);
    nind    = length(indd1);
    
    %no se encuentra todos los modos
    while nind<nmodej
        %no se borra los antiguos vectores, se guardan en DWN0.k2m1, detm1,
        %detfm1 (m1 por la iteracion corriente menos 1)
        DWN0.k2m1 	= DWN0.k2;
        detfm1      = detf;
        detm1       = det;
        
        %se inserta 1 punto entre cada punto del antiguo vector        
        DWN0.k2 	= 0.5*(DWN0.k2(1:end-1)+DWN0.k2(2:end));
        nk          = length(DWN0.k2);
        
        if  pol==1
            det = mode_Love(para,DWN0);
        else
            det = mode_Rayleigh_2(para,DWN0);
        end
        det(det==inf)=0;
        detf    = angle(det)-pi/2;
        
        %se reincorpora los nuevos vectores en los antiguos
        DWN0.k2m    = DWN0.k2;
        detfm       = detf;
        detm        = det;
        
        DWN0.k2(1:2:(2*nk+1)) 	= DWN0.k2m1;
        DWN0.k2(2:2:(2*nk)) 	= DWN0.k2m;
        detf(1:2:(2*nk+1))      = detfm1;
        detf(2:2:(2*nk))        = detfm;
        det(1:2:(2*nk+1))       = detm1;
        det(2:2:(2*nk))         = detm;
        

        nk          = 2*nk+1;
        indk        = 1:(nk-1);
        indd1       = logical(detf(1:(nk-1)).*detf(2:nk)<=0);
        indd1       = indk(indd1);
        nind        = length(indd1);
        
%         if nk/nk0>100
%             DWN0.omegac = wj+dw/100;
%         end
        if nk/nk0>200
            if nind==nmodej-1
                indd1       =[1,indd1];
                nind        = nind+1;
                det         =[-det(1);det];
            else
                pbw(j)      = 1;
                %se completa los vector con el numero de modos faltantes
                %para ensurar la convergencia al error propuesto para los
                %modos encontrados
                disp('pb faltan modos, hacer interpolacion')
                nmodej=nind;
%                 [k0f,indf,nmodej,det0]=solve_modos_falt_u2d(nf,j,wi,k01,DWN0,indd1,nmodej,para,ki);
                %reemplazar en el antiguo vector
%                 for i=1:length(indf)
            end
            break
        end
    end
    
    if  nind>nmodej
        if nind==nmodej+1 && sum((wj-w0c>=-1e-2*(w0c(2)-w0c(1))).*(wj-w0c<0))
            %problema de identificacion del primer punto del mayor modo
            %corriente # nmodej
            nmodej=nmodej+1;
        else
            disp('¿demasiados modos?')
            continue
        end
    end
    
    %mejora la convergencia pero estaria mucha mas sencillo en
    %interpolation lin
    det0  	= real(det(indd1));
    det1  	= real(det(indd1+1));
    k200  	= DWN0.k2(indd1);
    k201  	= DWN0.k2(indd1+1);
    
    k20up   = 0.5*(k200+k201).';
    vp0     = wj./k200;
    vp1     = wj./k201;
    
    %se recalcula puntos por los cuales el error entre el punto
    %encontrado y el punto interpolado es mayor al error
    while sum(abs(vp1-vp0)>err)>0
        k20     = k20up;
        DWN1    = DWN0;
        DWN1.k2	= k20.';
        if  pol==1
            det = mode_Love(para,DWN1);
        else
            det = mode_Rayleigh_2(para,DWN1);
        end
        det(det==inf)=0;
        
        %escoger el det de mismo signo y reemplazar
        ind01   = (sign(det0)==sign(real(det)));
        det0( ind01)=det(ind01);
        det1(~ind01)=det(~ind01);
        k200( ind01)=k20(ind01);
        k201(~ind01)=k20(~ind01);
        
        k20up   = 0.5*(k200+k201).';
        
        
        vp0     = wj./k200;
        vp1     = wj./k201;
        
    end
    a       = (det1-det0)./(k201-k200).';
    b       =  det1-a.*k201.';
    k20     = -b./a;
    
    %falto modos!
    if pbw(j) == 1
        k20         = [k20;zeros(nmodej-length(indd1),1)];
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
    plot(kj,wj,'c--')
    
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
    plot(f1(j,1:ikmax(j)),vp(j,1:ikmax(j)),'k--')%
end
xlabel('f');ylabel('vp');

figure(207);hold on;
for j=1:nmode
    plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'--')%
end
xlabel('f');ylabel('vg');