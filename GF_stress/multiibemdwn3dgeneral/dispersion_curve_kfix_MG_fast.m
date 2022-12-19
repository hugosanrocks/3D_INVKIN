function [k2,w0,ikmax]=dispersion_curve_kfix_MG_fast(para)

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


tic

ncapas  = para.nsubmed;
para    = Vp2Cij(para);
for ms=1:ncapas
    para.reg(1).sub(ms).Ci 	= para.reg(1).sub(ms).C;
    para.reg(1).sub(ms).tipoatts 	= 3;
end
wmax    = 2*pi*para.fmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour k=0, il n y a que des modes de volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vT      = zeros(ncapas,1);
for im  = 1:ncapas-1
    C       = para.reg(1).sub(im).C;
    rho     = para.reg(1).sub(im).rho;
    vT(im) 	= sqrt(C(6,6)/rho);
end
im      = ncapas;
C       = para.reg(1).sub(im).C;
rho     = para.reg(1).sub(im).rho;
vT(im) 	= sqrt(C(6,6)/rho);
vdeep   = vT(end)*(1-1e-6);

vTtmp           = vT;
vTtmp(vT>=vdeep)= [];
vdeep2          = max(vTtmp);

% if ncapas == 2
%     imode=1;
%     w00=0;
%     k20=0;
%     while w00(imode)<wmax
%         k20(imode+1)  = vT(1)*imode*pi/(h(1)*sqrt(vT(2)^2-vT(1)^2));
%         w00(imode+1)  = k20(imode+1)*vT(2);
%         imode       = imode+1;
%     end
%     k20=k20(1:(imode-1));
%     w00=w00(1:(imode-1));
%     nmode=imode-1;
% end

[k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,vdeep);
nmode       = length(k20);

figure(205);hold on;plot(k20,w00,'.r')

k2      = zeros(1000,nmode);
w0      = zeros(1000,nmode);
err1    = zeros(1000,nmode);
ikmax 	= zeros(1,nmode);
w0(1,:)	= w00;
k2(1,:)	= k20;

facDK   = 20;
for imode=1:nmode
    disp(imode);
    ik  = 1;
    if imode==1 && nmode==1
        DK = (k2(1,imode))/facDK;
    elseif imode+1<=nmode
        DK = (k2(1,imode+1)-k2(1,imode))/facDK;
    else
        DK = (k2(1,imode)-k2(1,imode-1))/facDK;
    end
    DK1     = facDK*DK;
    %     while k2(ik,imode)<kmax
    while w0(ik,imode)<wmax
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
            nextw0  = k2(ik1,imode)*vdeep;
        end
        
        %premiere approx de l intervale de recherche
        if ik1>4
            dw0=2*err*DK;
            dw0(dw0==0)=.01;
        else
            dw0 = abs(w0(ik,imode)-nextw0)/10;
            dwmin=k2(ik1,imode)*(vdeep-vdeep2)/2;
            dw0 =min(dw0,dwmin);
            %pourt le premier point
            if nextw0-dw0<0
                dw0=nextw0;
            end
        end
        
        %verifier que la partie reelle et imaginaire change de signe
        wi  = nextw0-dw0;
        wf  = nextw0+dw0;
        [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        %         plot(k2(ik1,imode),nextw0,'c.')
        
        %verifier qu on n est pas proche d une branche asymptotique
        %ds ce cas la il y a une des parties qui croise 0 mais rechange
        %vite de signe, l intervalle de recherche doit etre tres petit
        if (sign0>0 && sign1<=0 || sign0<=0 && sign1>0) && ik1>4
            if min(abs(nextw0./vT-k2(ik1,imode)))<DK
                %                 plotthepb(wi,wf,para,DWN)
                
                if sign0>0 && sign1<=0
                    %chercher le croisement de imag
                    w00   = fzero(@(w) zerodeAwIR(para,DWN,w,2),[wi wf]);
                    %verifier le chgmt de sign de l autre partie
                    [sign0,~]=checkchgmtsign(w00-dw0*1e-3,w00+dw0*1e-3,para,DWN);
                    if sign0<=0
                        wi  = w00-dw0*1e-3;
                        wf  = w00+dw0*1e-3;
                    end
                else
                    %chercher le croisement de imag
                    w00   = fzero(@(w) zerodeAwIR(para,DWN,w,1),[wi wf]);
                    %verifier le chgmt de sign de l autre partie
                    [~,sign1]=checkchgmtsign(w00-dw0*1e-3,w00+dw0*1e-3,para,DWN);
                    if sign1<=0
                        wi  = w00-dw0*1e-3;
                        wf  = w00+dw0*1e-3;
                    end
                end
            end
        end
        
        dw00=dw0;
        nextw00=nextw0;
        %d abord tenter de reduire avant d agrandir
        %pb qd l intervale de recherche est trop grand ou qu on voit 2 zeros
        j=1;
        while (sign0>=0 || sign1>=0) && j<6
            % plotthepb(wi,wf,para,DWN)
            %             if sign0<0
            dw=wf-wi;
            [sign0]=checkchgmtsign(wi,wi+dw/2,para,DWN);
            if sign0<0
                nextw0=wi+dw/4;
            else
                nextw0=wi+3*dw/4;
            end
            %             else
            %                 dw0 = dw0*4;
            %             end
            j   = j+1;
            dw0 = dw0/2;
            wi  = nextw0-dw0;
            wf  = nextw0+dw0;
            [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        end
        
        %pb qd l intervale de recherche est trop petit
        if sign0>=0 && sign1>=0
            % plotthepb(wi,wf,para,DWN)
            dw0     = dw00*1.5;
            nextw0  = nextw00;
            wi      = nextw0-dw0;
            wf      = nextw0+dw0;
            [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        end
        
        
        if sign0<=0 && sign1<=0
            %tt est ok on recherche
            w0(ik1,imode)   = fzero(@(w) zerodeAw(para,DWN,w),[wi wf]);
            
            %on verifie si il n y a pas un faux zero
            [sign0,sign1]=checkchgmtsign(w0(ik1,imode)-dw0*1e-3,w0(ik1,imode)+dw0*1e-3,para,DWN);
            if sign0<=0 && sign1<=0
                %                 if 2*min(abs(w0(ik1,imode)./vT-k2(ik1,imode)))<DK && ik>4
                %                     DK=DK/4;
                %                 end
                err             = abs(w0(ik1,imode)-nextw0)/DK;
                err1(ik1,imode) = err*DK;
                ik              = ik+1;
                DK              = min(1.3*DK,DK1);
                
                %                                 plot(k2(ik1,imode),w0(ik1,imode),'r.')
                %                                 ik1
            else
                DK              = DK/4;
                if ik1>4
                    err=err*2;
                end
            end
            %             toto=0;
        else
            %si rien n a marche avant, on diminue le pas d avance
            DK              = DK/2;
        end
%         plot(k2(ik,imode),w0(ik,imode),'c.')
        
    end
    
    ikmax(imode)=ik1;
    %     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'c')
    %     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')
    
end
toc

% figure(205);hold on;
% for imode=1:nmode
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
%     plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')
% end

% figure(1);hold on;
% for imode=1:nmode
%     plot3(k2(1:ikmax(imode),imode)/126*620,w0(1:ikmax(imode),imode)/126*401,20+w0(1:ikmax(imode),imode)*0,'r')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vp=w/k =f(w)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(206);hold on;
for imode=1:nmode
    vp=w0(1:ikmax(imode),imode)./k2(1:ikmax(imode),imode);
    if imode==1
        vp(1)=vdeep;
    end
    plot(w0(1:ikmax(imode),imode)/2/pi,vp,'k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(207);hold on;
for imode=1:nmode
    vg=diff(w0(1:ikmax(imode),imode))./diff(k2(1:ikmax(imode),imode));
    vg=[vdeep;vg]; %#ok<AGROW>
    f1=[w0(1,imode);(w0(1:ikmax(imode)-1,imode)+w0(2:ikmax(imode),imode))/2]/2/pi;
    plot(f1,vg,'k')
end