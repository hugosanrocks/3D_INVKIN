function [vg,f1,k2,w0,ikmax]=dispersion_curve_kfix_MG_fast2(para)

% cambios con respecto a  function dispersion_curve_kfix(para):
% en lugar de buscar los zeros por cada k en un gran intervalo,
% se aprovecha el conocimiento de los zeros anteriores para reducir los
% intervalos de busqueda siguiendo cada modo
% los primeros zeros se calculan con asintotas en el caso 1 estrato /SE o
% con busqueda numerica sobre la asintota w=k.cmax (amejorar)
% pb: se puede equivocar con otros modos muy cercano y pasarse de uno a
% otro o equivocarse con las asintotas que son tambien solucion (sin cambio
% de signo pero con un 0)
% cambios con respecto a  function dispersion_curve_kfix_MG_fast(para):
% l intervale de recherche est amélioré en prenant en compte l avance du DK
% et l erreur a la tentative anterieure
% la recherche ne se fait plus avec fzero pour aller plus vite
% on ralentit le pas Dk devant les branches asymptotique, cette fois ci les
% branches asymptotiques sont correctement détectées
% le gain de tps est a peu pres d un factuer 2

ncapas  = para.nsubmed;
para    = Vp2Cij(para);

for ms=1:ncapas
    para.reg(1).sub(ms).Ci 	= para.reg(1).sub(ms).C;
end

wmax    = 2*pi*para.fmax;
facav   = 2;


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
vmax    = max(vT);
vTtmp   = vT;
vTtmp(vT==vmax)=[];
vmax2   = (1-1e-8)*max(vTtmp);
vmax    = (1-1e-8)*vmax;

tic

[k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,vmax);
nmode       = length(k20);
k2          = zeros(1000,nmode);
w0          = zeros(1000,nmode);
err1        = zeros(1000,nmode);
ikmax       = zeros(1,nmode);
w0(1,:)     = w00;
k2(1,:)     = k20;


facDK   = 20;
DK0     = 0;
for imode=1:nmode
    if imode==1 && nmode==1
        DK = .1;
    elseif imode+1<=nmode
        DK = (k2(1,imode+1)-k2(1,imode));
    else
        DK = (k2(1,imode)-k2(1,imode-1));
    end
    DK0=DK0+DK;
end
DK0= DK0/nmode/facDK;
DK1= facDK*DK0;

for imode=1:nmode
%     disp(imode);
    facDK   = 20;
    ik      = 1;
    DK      = DK0;
    w0p     = 0;
    
    %     while k2(ik,imode)<kmax
    while (w0(ik,imode)<wmax) || w0p==0
        if w0(ik,imode)>wmax
            %pour avoir wmax
            w0p=1;
        end
        ik1     = ik+1;
        
        %avance de k2
        k2(ik1,imode)	= k2(ik,imode) + DK;
        %         if k2(ik1,imode)>kmax
        %             k2(ik1,imode)=kmax;
        %         end
        DWN.k2      = k2(ik1,imode);
        
        %approx de la solution par extrapolation
        if ik>=4
            nextw0 = interp1(k2(1:ik,imode),w0(1:ik,imode),k2(ik1,imode),'pchip','extrap');
        elseif ik==3
            w1      = w0(1:3,imode);
            k21     = k2(1:3,imode);
            dwdk    = ((w1(2)-w1(3))*(k21(1)-k21(3))^2-(w1(1)-w1(3))*(k21(2)-k21(3))^2)/...
                ((k21(2)-k21(3))*(k21(1)-k21(3))^2-(k21(1)-k21(3))*(k21(2)-k21(3))^2);
            d2wdk2  = ((w1(2)-w1(3))-dwdk*(k21(2)-k21(3)))*2/(k21(2)-k21(3))^2;
            nextw0  = w1(3)+dwdk*(k2(ik1,imode)-k2(ik,imode))+d2wdk2*(k2(ik1,imode)-k2(ik,imode))^2/2;
        elseif ik==2
            dwdk    = (w0(1,imode)-w0(2,imode))/(k2(1,imode)-k2(2,imode));
            nextw0  = w0(2,imode)+dwdk*(k2(ik1,imode)-k2(ik,imode));
        elseif ik==1
            nextw0  = k2(ik1,imode)*vmax;
        end
        
        %premiere approx de l intervale de recherche
        if ik1>4
            dw0=err1(ik,imode)*facDK*DK;
            dw0(dw0<1e-6)=1e-6;%pb precision machine
        else
            dw0 = abs(w0(ik,imode)-nextw0)/2;
            dwmin=k2(ik1,imode)*(vmax-vmax2)/2;
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
        %         plot(k2(ik1,imode),nextw0,'c.')
        
        dw00=dw0;
        
        %verifier qu on n est pas proche d une branche asymptotique
        %ds ce cas la il y a une des parties qui croise 0 mais rechange
        %vite de signe, l intervalle de recherche doit etre tres petit
        %qd le sign==0 alors c est que c est des ondes homogenes
        if (sign0>0 && sign1<=0 || sign0<=0 && sign1>0) || (min(abs(nextw0./vT-k2(ik1,imode)))<2*DK)
            %                 plotthepb(wi-dw0,wf+dw0,para,DWN)
            npt     = min(max(fix(abs(dw0)*2/1e-4),10),50);
            [~,~,det1,para1,DWN1]=checkchgmtsign_vec(wi-dw0,wf+dw0,para,DWN,npt);
            %   figure(203);hold on;plot(DWN1.omegac,real(det1),'k');plot(DWN1.omegac,imag(det1),'r');plot(DWN1.omegac,abs(det1),'')
            np=min(max(fix(abs(dw0/npt)*2/1e-6),10),100);
            w0000=cherchermin_k(abs(det1),para1,DWN1,np);
            nind=length(w0000);
            if nind~=0
                %on elimine les asymptotes
                dwl=(DWN1.omegac(3)-DWN1.omegac(1))/np*2;
                if para.pol==2
                    while nind~=0
                        if (min(min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vP.'*k2(ik1,imode))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vP.'*k2(ik1,imode)),[],2);
                            [~,ind0]=min(abs(w0000-vP(ind0)*k2(ik1,imode)));
                            w0000(ind0)=[];
                            nind=length(w0000);
                        else
                            break
                        end
                    end
                end
                if nind~=0
                    while nind~=0
                        if (min(min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vT.'*k2(ik1,imode))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(w0000*ones(1,ncapas)-ones(nind,1)*vT.'*k2(ik1,imode)),[],2);
                            [~,ind0]=min(abs(w0000-vT(ind0)*k2(ik1,imode)));
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
                    w0(ik1,imode)   = w0000;
                    err             = abs(w0(ik1,imode)-nextw0)*2/DK;
                    err             = max(err,err1(ik,imode)/2);
                    err1(ik1,imode) = err;
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
            npt=min(max(fix(abs(dw0)*2/1e-4),10),50);
            [~,~,det1,para1,DWN1]=checkchgmtsign_vec(wi,wf,para,DWN,npt);
            np=min(max(fix(abs(dw0/npt)*2/1e-6),10),100);
            w0000=cherchermin_k(abs(det1),para1,DWN1,np);
            nind=length(w0000);
            if nind~=1
                DK              = DK/2;
                continue
            else
                w0(ik1,imode)   = w0000;
                if j==1
                    err             = abs(w0(ik1,imode)-nextw0)*2/DK;
                    err             = max(err,err1(ik,imode)/2);
                    err1(ik1,imode) = err;
                    DK              = min(facav*DK,DK1);
                else
                    err             = dw0/DK;
                    err1(ik1,imode) = err;
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
    % en caso de problema
    % plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
    % plot(k2(1:(ik1),imode),[w0(1:(ik),imode);nextw0],'r.')
    % plot(k2(1:ik,imode),w0(1:ik,imode)+ 2*err1(1:ik,imode),'r')
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%acpopt
figure(205);hold on;
plot(k20,w00,'.r')
for imode=1:nmode
    %     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
    %     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')
end
xlabel('k_x');
ylabel('\omega');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica asintotas                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%acpopt
k2as  = 0:1e3;
for ic=1:ncapas
    C   = para.reg(1).sub(ic).C;
    rho = para.reg(1).sub(ic).rho;
    vT  = sqrt(C(6,6)/rho);
    vL  = sqrt(C(1,1)/rho);
    if  para.pol==1
        w=k2as*vT;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2as(1:i0),w(1:i0),'k');
    else
        w=k2as*vT;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2as(1:i0),w(1:i0),'');
        w   = k2as*vL;
        i0=find(w>wmax,1);
        figure(205);hold on;plot(k2as(1:i0),w(1:i0),'r');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vp=w/k =f(w)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%acpopt
figure(206);hold on;
for imode=1:nmode
    vp=w0(1:ikmax(imode),imode)./k2(1:ikmax(imode),imode);
    if imode==1
        vp(1)=vmax;
    end
    plot(w0(1:ikmax(imode),imode)/2/pi,vp,'k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(207);hold on;%acpopt (a commenter pour l optimisation)
vg=zeros(nmode,max(ikmax)-1);
f1=zeros(nmode,max(ikmax)-1);

for imode=1:nmode
    vg(imode,1:(ikmax(imode)-1))=diff(w0(1:ikmax(imode),imode))./diff(k2(1:ikmax(imode),imode));
    f1(imode,1:(ikmax(imode)-1))=(w0(1:ikmax(imode)-1,imode)+w0(2:ikmax(imode),imode))/2/2/pi;
    plot(f1(imode,1:(ikmax(imode)-1)),vg(imode,1:(ikmax(imode)-1)),'k')%acpopt
end
ikmax=ikmax-1;