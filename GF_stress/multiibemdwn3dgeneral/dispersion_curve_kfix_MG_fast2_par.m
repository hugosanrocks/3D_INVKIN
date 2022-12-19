function [vg,f1,ikmax]=dispersion_curve_kfix_MG_fast2_par(para)
% OJO !!!!! no actualizado con respecto a dispersion_curve_kfix_MG_fast2
% pero la philisofia esta aqui
% cambios con respecto a  function dispersion_curve_kfix(para):
% en lugar de buscar los zeros por cada k en un gran intervalo,
% se aprovecha el conocimiento de los zeros anteriores para reducir los
% intervalos de busqueda siguiendo cada modo
% los primeros zeros se calculan con asintotas en el caso 1 estrato /SE o
% con busqueda numerica sobre la asintota w=k.cmax (amejorar)
% pb: se puede equivocar con otros modos muy cercano y pasarse de uno a
% otro o equivocarse con las asintotas que son tambien solucion (sin cambio
% de signo pero con un 0)
% se corrige este problema en la version xxx donde se precisa el intervalo
% de busqueda procesando varios modos al mismo k


wmax    = 2*pi*para.fmax;
facav   = 2;
ncapas  = para.nsubmed;
para    = Vp2Cij(para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour k=0, il n y a que des modes de volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h       = zeros(ncapas,1);
vT      = zeros(ncapas,1);
vP      = zeros(ncapas,1);

for im  = 1:ncapas-1
    h(im)   = para.reg(1).sub(im).h;
    C       = para.reg(1).sub(im).C;
    rho     = para.reg(1).sub(im).rho;
    vT(im) 	= sqrt(C(6,6)/rho);
    vP(im) 	= sqrt(C(1,1)/rho);
end
im      = ncapas;
C       = para.reg(1).sub(im).C;
rho     = para.reg(1).sub(im).rho;
vT(im) 	= sqrt(C(6,6)/rho);
vP(im) 	= sqrt(C(1,1)/rho);
vmax    =max(vT);
vTtmp   =vT;
vTtmp(vT==vmax)=[];
vmax2   =max(vTtmp);


if ncapas ==2 && para.pol==1
    imode=1;
    w=0;
    while w<wmax
        k20(imode+1)= vT(1)*imode*pi/(h(1)*sqrt(vT(2)^2-vT(1)^2));
        w00(imode+1)= k20(imode+1)*vT(2);
        w         	= w00(imode+1);
        imode       = imode+1;
        %     figure(205);hold on;plot(k2(1,imode+1),w0(1,imode+1),'.k')
    end
    k20=k20(1:(imode-1));
    w00=w00(1:(imode-1));
    nmode=imode-1;
else
    if  para.pol==2
        vmax=(1-1e-8)*vmax;
    end
    [k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,vmax);
    nmode       = length(k20);
end

k2      = zeros(nmode,1000);
w0      = zeros(nmode,1000);
err1    = zeros(nmode,1000);
ikmax 	= zeros(nmode,1);
w0(:,1)	= w00;
k2(:,1)	= k20;

k2par   = k2(1:nmode,1);%for parallelization
w0par   = w0(1:nmode,1);

facDK0   = 20;
DK0     = 0;
for imode=1:nmode
    if imode==1 && nmode==1
        DK = .1;
    elseif imode+1<=nmode
        DK = (k2(imode+1,1)-k2(imode,1));
    else
        DK = (k2(imode,1)-k2(imode-1,1));
    end
    DK0=DK0+DK;
end
DK0= DK0/nmode/facDK0;
DK1= facDK0*DK0;

parfor imode=1:nmode
    facDK   = facDK0;
    errp	= zeros(1000,1);%for parallelization
    w0p     = zeros(1000,1);%for parallelization
    k2p     = zeros(1000,1);%for parallelization
    k2p(1)  = k2par(imode);
    w0p(1)  = w0par(imode);
    
    ik      = 1;
    ik1     = 0;
    nextw0  = 0;
    DK      = DK0;
    w0b     = 0;
    %     while k2(ik)<kmax
    while (w0p(ik)<wmax) || w0b==0
        if w0p(ik)>wmax
            %pour avoir wmax
            w0b=1;
        end
        ik1     = ik+1;
        
        %avance de k2
        k2p(ik1)	= k2p(ik) + DK;
        %         if k2p(ik1)>kmax
        %             k2p(ik1)=kmax;
        %         end
        DWN = struct('k2',k2p(ik1));

        
        %approx de la solution par extrapolation
        if ik>=4
            nextw0 = interp1(k2p(1:ik),w0p(1:ik),k2p(ik1),'pchip','extrap');
        elseif ik==3
            w1      = w0p(1:3);
            k21     = k2p(1:3);
            dwdk    = ((w1(2)-w1(3))*(k21(1)-k21(3))^2-(w1(1)-w1(3))*(k21(2)-k21(3))^2)/...
                ((k21(2)-k21(3))*(k21(1)-k21(3))^2-(k21(1)-k21(3))*(k21(2)-k21(3))^2);
            d2wdk2  = ((w1(2)-w1(3))-dwdk*(k21(2)-k21(3)))*2/(k21(2)-k21(3))^2;
            nextw0  = w1(3)+dwdk*(k2p(ik1)-k2p(ik))+d2wdk2*(k2p(ik1)-k2p(ik))^2/2;
        elseif ik==2
            dwdk    = (w0p(1)-w0p(2))/(k2p(1)-k2p(2));
            nextw0  = w0p(2)+dwdk*(k2p(ik1)-k2p(ik));
        elseif ik==1
%             if imode==1
                nextw0  = k2p(ik1)*vmax;
%             else
%                 nextw0  = w0p(1)+ (w0p(imode)-w0p(1))/(k2p(imode)-k2p(1))*(k2p(ik1)-k2p(1));
%             end
        end
        
        %premiere approx de l intervale de recherche
        if ik1>4
            dw0=errp(ik)*facDK*DK;
            dw0(dw0<1e-6)=1e-6;%pb precision machine
        else
            dw0 = abs(w0p(ik)-nextw0)/2;
            dwmin=k2p(ik1)*(vmax-vmax2)/2;
            dw0 =min(dw0,dwmin);
            %pourt le premier point
            if nextw0-dw0<0
                dw0=nextw0;
            end
        end
        
        %verifier que la partie reelle et imaginaire change de signe
        wi  = nextw0-dw0;
        wf  = nextw0+dw0;
        [sign0,sign1]=checkchgmtsign_vec(wi,wf,para,DWN,2);
        %         plot(k2p(ik1),nextw0,'c.')
        
        dw00=dw0;
        
        %verifier qu on n est pas proche d une branche asymptotique
        %ds ce cas la il y a une des parties qui croise 0 mais rechange
        %vite de signe, l intervalle de recherche doit etre tres petit
        if (sign0>0 && sign1<=0 || sign0<=0 && sign1>0) || (min(abs(nextw0./vT-k2p(ik1)))<2*DK)
            %                 plotthepb(wi-4*dw0,wf+4*dw0,para,DWN)
            
            npt=min(max(fix(abs(dw0)*2/1e-4),3),50);
            [~,~,det1,para1,DWN1]=checkchgmtsign_vec(wi-dw0,wf+dw0,para,DWN,npt);
            %   figure(203);hold on;plot(DWN1.omegac,real(det1),'k');plot(DWN1.omegac,imag(det1),'r');plot(DWN1.omegac,abs(det1),'')
            w0000=cherchermin_k(abs(det1),para1,DWN1,100);
            nind=length(w0000);
            if nind~=0
                %on elimine les asymptotes
                dwl=(DWN1.omegac(3)-DWN1.omegac(1))/100*2;
                if para.pol==2
                    while nind~=0
                        if (min(min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vP.'*k2p(ik1))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vP.'*k2p(ik1)),[],2);
                            [~,ind0]=min(abs(w0000-vP(ind0)*k2p(ik1)));
                            w0000(ind0)=[];
                            nind=length(w0000);
                        else
                            break
                        end
                    end
                end
                if nind~=0
                    while nind~=0
                        if (min(min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vT.'*k2p(ik1))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vT.'*k2p(ik1)),[],2);
                            [~,ind0]=min(abs(w0000-vT(ind0)*k2p(ik1)));
                            w0000(ind0)=[];
                            nind=length(w0000);
                        else
                            break
                        end
                    end
                end
                if nind>1
                    DK=DK/2;
                    if ik>5
                        ik = ik-1;
                    end
                    continue
                elseif nind~=0
                    w0p(ik1)   = w0000;
                    err             = abs(w0p(ik1)-nextw0)*2/DK;
                    err             = max(err,errp(ik)/2);
                    errp(ik1) = err;
                    ik              = ik+1;
                    continue
                end
                facDK   = 1;
            end
        end
        
        %pb qd l intervale de recherche est trop petit
        j=1;
        dw0=dw00;
        while (sign0>0 || sign1>0) && j<4
            % plotthepb(wi,wf,para,DWN)
            dw0 = dw0*2;
            wi  = nextw0-dw0;
            wf  = nextw0+dw0;
            [sign0,sign1,~,~,~]=checkchgmtsign_vec(wi,wf,para,DWN,2);
            j=j+1;
        end
        
        if sign0<=0 && sign1<=0
            %tt est ok on recherche
            npt=min(max(fix(abs(dw0)*2/1e-4),3),50);
            [~,~,det1,para1,DWN1]=checkchgmtsign_vec(wi,wf,para,DWN,npt);
            w0000=cherchermin_k(abs(det1),para1,DWN1,100);
            nind=length(w0000);
            if nind~=1
                DK              = DK/2;
                continue
            else
                w0p(ik1)   = w0000;
                if j==1
                    err             = abs(w0p(ik1)-nextw0)*2/DK;
                    err             = max(err,errp(ik)/2);
                    errp(ik1) = err;
                    DK              = min(facav*DK,DK1);
                else
                    err             = dw0/DK;
                    errp(ik1) = err;
                end
                ik              = ik+1;
            end
            facDK   = 1;
        else
            %si rien n a marche avant, on diminue le pas d avance
            DK      = DK/2;
            facDK   = j;
        end
        
    end
    ikmax(imode)=ik1;
    k2(imode,:) = k2p;
    w0(imode,:) = w0p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vg=zeros(nmode,max(ikmax)-1);
f1=zeros(nmode,max(ikmax)-1);

for imode=1:nmode
    vg(imode,1:(ikmax(imode)-1))=diff(w0(imode,1:ikmax(imode)))./diff(k2(imode,1:ikmax(imode)));
    f1(imode,1:(ikmax(imode)-1))=(w0(imode,1:ikmax(imode)-1)+w0(imode,2:ikmax(imode)))/2/2/pi;
end
ikmax=ikmax-1;