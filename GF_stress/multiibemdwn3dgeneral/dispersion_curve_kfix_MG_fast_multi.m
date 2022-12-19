function [k2,w0,ikmax]=dispersion_curve_kfix_MG_fast_multi(para)

% cambios con respecto a "function dispersion_curve_kfix_MG_fast(para)":
% en lugar de calcular por cada modo en funcion de k, 
% se calcula todos los modos en function de k hasta la frecuencia max,al fin de
% aprovechar el conocimiento de los otros modos para mejor discriminar el
% intervalo de busqueda de los zeros
% no mejora el tiempo pero si el intervalo de busqueda

% cambios con respecto a  function dispersion_curve_MG_kfix(para):
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
end
wmax    = 2*pi*para.fmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour k=0, il n y a que des modes de volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncapas  = para.nsubmed;
h       = zeros(ncapas,1);
vT      = zeros(ncapas,1);

for im  = 1:ncapas-1
    h(im)   = para.reg(1).sub(im).h;
    C       = para.reg(1).sub(im).C;
    rho     = para.reg(1).sub(im).rho;
    vT(im) 	= sqrt(C(6,6)/rho);
end
im      = ncapas;
C       = para.reg(1).sub(im).C;
rho     = para.reg(1).sub(im).rho;
vT(im) 	= sqrt(C(6,6)/rho);
vmax    =max(vT);
vTtmp   =vT;
vTtmp(vT==vmax)=[];
vmax2   =max(vTtmp);


[k20,w00]   = dispersion_curve_k_critik_MG(para,wmax,vmax);
nmode       = length(k20);
figure(205);hold on;plot(k20,w00,'.r')

k2      = zeros(1000,nmode);
w0      = zeros(1000,nmode);
err1    = zeros(1000,nmode);
ikmax 	= zeros(1,nmode);
w0(1,:)	= w00;
k2(1,:)	= k20;
nmk2    = size(k2);
facDK   = 20;
ik      = ones(nmode,1);
nmodek  = 1;
imode   = 1:nmodek;
nextw0	= zeros(nmodek,1);


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

while w0(ik(1),1)<wmax
    % se toma en cuenta otro modo
    % y avance de todos los k2
    if k2(ik(1),1)+DK > k2(1,nmodek+1)
        nmodek          = nmodek+1;
        DK              = min(DK1/facDK,DK);
        imode           = 1:nmodek;
        ik1             = ik(imode) + 1;
        linearInd       = sub2ind(nmk2, ik1,imode.');
        k2(linearInd)	= k2(1,nmodek)+DK;
        nextw0          = zeros(nmodek,1);
    else
        ik1            	= ik(imode) + 1;
        linearInd       = sub2ind(nmk2, ik1,imode.');
        k2(linearInd)	= k2(ik(1),1)+DK;
    end

    if imode+1<nmode
        DK = (k2(1,imode+1)-k2(1,imode))/facDK;
    else
        DK = (k2(1,imode)-k2(1,imode-1))/facDK;
    end
    DK1     = facDK*DK*30;
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
            nextw0  = k2(ik1,imode)*vmax;
        end
        
        %premiere approx de l intervale de recherche
        if ik1>4
            dw0=2*err;
            dw0(dw0==0)=.01;
        else
            dw0 = abs(w0(ik,imode)-nextw0)/10;
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
        [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
%         plot(k2(ik1,imode),nextw0,'c.')
        
        %pb qd l intervale de recherche est trop petit
        j=1;
        while sign0>=0 && sign1>=0 && j<6
            % plotthepb(wi,wf,para,DWN)
            
            j   = j+1;
            dw0 = dw0*2;
            wi  = nextw0-dw0;
            wf  = nextw0+dw0;
            [sign0,sign1]=checkchgmtsign(wi,wf,para,DWN);
        end
        
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
        
        if sign0<=0 && sign1<=0
            %tt est ok on recherche
            ik  = ik+1;
            w0(ik1,imode)   = fzero(@(w) zerodeAw(para,DWN,w),[wi wf]);
            DK              = min(1.3*DK,DK1);
            err             = abs(w0(ik1,imode)-nextw0);
            err1(ik1,imode) = err;
%                      plot(k2(ik1,imode),w0(ik1,imode),'r.')
        else
            %si rien n a marche avant, on diminue le pas d avance
            DK              = DK/2;
        end
    end
    ikmax(imode)=ik1;
        plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')

end
toc

figure(205);hold on;
for imode=1:nmode
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)- 2*err1(1:ikmax(imode),imode),'r')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode)+ 2*err1(1:ikmax(imode),imode),'r')
end

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
        vp(1)=vmax;
    end
    plot(w0(1:ikmax(imode),imode)/2/pi,vp,'k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(207);hold on;
for imode=1:nmode
    vg=diff(w0(1:ikmax(imode),imode))./diff(k2(1:ikmax(imode),imode));
    vg=[vmax;vg]; %#ok<AGROW>
    f1=[w0(1,imode);(w0(1:ikmax(imode)-1,imode)+w0(2:ikmax(imode),imode))/2]/2/pi;
    plot(f1,vg,'k')
end