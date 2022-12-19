function [vg,f1,ikmax]=dispersion_curve_kfix_Haskel_fast_check_saut(para)

% cambios con respecto a  function dispersion_curve_kfix(para):
% en lugar de buscar los zeros por cada k en un gran intervalo,
% se aprovecha el conocimiento de los zeros anteriores para reducir los
% intervalos de busqueda siguiendo cada modo
% los primeros zeros se calculan con asintotas en el caso 1 estrato /SE o
% con busqueda numerica sobre la asintota w=k.cmax (amejorar)
% pb: se puede equivocar con otros modos muy cercano y pasarse de uno a
% otro o equivocarse con las asintotas que son tambien solucion (sin cambio
% de signo pero con un 0)
% Se puede mejorar esta version con un mex pero tomaria mucho tiempo.

para.DWNomei=0;

wmax    = 2*pi*para.fmax;

facav   = 1.3;
ncapas  = para.nsubmed;
para    = Vp2Cij(para);
para.sub= para.reg(1).sub;

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
vmax2   = (1-1e-8)*max(vTtmp);
vmax    = (1-1e-8)*vmax;
kmax    = wmax/min(vT);

[k20,w00]   = dispersion_curve_k_critik_Haskel(para,wmax,vmax);
nmode       = length(k20);

vasin   = unique(vT);
nm      = length(vasin);
k2as    = zeros(nmode,nm);
w0as    = zeros(nmode,nm);
k2as(:,nm)=k20;
w0as(:,nm)=w00;
for iasin=2:(nm-1)
    [ktmp,wtmp]=dispersion_curve_k_critik_Haskel(para,2*wmax,vasin(iasin));
    ni=min(length(ktmp),nmode);
    k2as(1:ni,iasin)=ktmp(1:ni);
    w0as(1:ni,iasin)=wtmp(1:ni);
end


figure(205); hold on

k2      = zeros(1000,nmode);
w0      = zeros(1000,nmode);
err1    = zeros(1000,nmode);
ikmax 	= zeros(1,nmode);
w0(1,:)	= w00;
k2(1,:)	= k20;


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
DK0= DK0/nmode/10;
DK1= 100*DK0;%paso max

tic
% nmode=1;
for imode=1:nmode
    imode
    facDw   = 1;%ouverture en Dw
    ik      = 1;
    DK      = DK0;
    
    while (w0(ik,imode)<wmax)
        %etape de verification pour eviter de passer d un mode a un autre
        %on regarde entre les 2 dernieres itérations pour voir si il y a un
        %zero a l endroit median
        if ik>5 && imode>1
            k2t     = mean(k2((ik-2):(ik-1),imode));
            DWN.k2  = k2t;
            w0c     = interp1(k2(ik-3:ik,imode),w0(ik-3:ik,imode),k2t,'pchip');
            w0cl(1) = interp1(k2(ik-2:ik-1,imode),w0(ik-2:ik-1,imode),k2t,'linear');
            w0cl(2) = interp1(k2(ik-2:2:ik,imode),w0(ik-2:2:ik,imode),k2t,'linear');
            w0cl(3) = interp1(k2(ik-3:2:ik-1,imode),w0(ik-3:2:ik-1,imode),k2t,'linear');
            %             w0cl(4) = interp1(k2(ik-3:3:ik,imode),w0(ik-3:3:ik,imode),k2t,'linear');
            dw0     = max(abs(w0c-w0cl(1)));
            wi      = w0c-dw0;
            wf      = w0c+dw0;
            
            %             [~,~,~,det1,para1,DWN1]=checkchgmtsign_H_vec(wi,wf,para,DWN,2);
            %             solw    = cherche_zero_dic(para1,DWN1,det1);
            %             nind    = length(solw);
            %             if nind~=1
            %plotthepb_H(wi,wf,para,DWN)
            
            signp   = checkchgmtsign_H_vec(wi,wf,para,DWN,2);
            if signp>0
                dw0     = max(abs(w0c-w0cl));
                wi      = w0c-dw0;
                wf      = w0c+dw0;
                signp   = checkchgmtsign_H_vec(wi,wf,para,DWN,2);
                if signp>0
                    ik  = ik-2;
                    DK  = DK/2;
                end
            end
        end
        
        ik1     = ik+1;
        
        %avance de k2
        %on verifie la proximite d une branche asymptotique qd il faut la
        %croiser, on exlue donc la derniere branche asymptotique (1 ere
        %couche)
        [~,indmv]   = min(abs(w0(ik,imode)-k2(ik,imode).*vT(2:end)));
        dpb         = k2(ik,imode)*vT(1+indmv)-w0(ik,imode);
        if abs(dpb)<2*DK && dpb<0
            DK = min(DK1/10,DK);
        end
        %on verifie aussi la proximite des deviations rapides dues a des
        %croisements des branches asymptotiques qd les croisements sont
        %plus rapproches que l espace entre les points de departs
        %
        %identification du secteur indmv
        indmv   = find(w0(ik,imode)>k2(ik,imode).*vT(2:end),1,'last');
        if ~isempty(indmv)
            %le prochain embranchement du secteur
            na=find(w0(ik,imode)<w0as(:,indmv+1),1,'first');
            %si l espace au prochain branchement est plus proche que l
            %eloignement entre les 2 prochaines branches
            if ~isempty(na)
                tmp=w0as(na+1,indmv+1)-w0as(na,indmv+1);
                if tmp<0
                    tmp=w0as(na,indmv+1)-w0as(na-1,indmv+1);
                else
                    tmp=min(tmp,w0as(na,indmv+1)-w0as(na-1,indmv+1));
                end
                if abs(w0(ik,imode)-w0as(na,indmv+1))<tmp
                    %et si la distance des 2 prochains branchements est plus
                    %petites que 2 fois l espace entre les points de departs, il va
                    %y avoir des deviations rapides
                    if tmp<2*(w0as(2,nm)-w0as(1,nm))
                        DK = min(DK0,DK);
                    end
                end
            end
        end
        
        
        k2(ik1,imode)	= k2(ik,imode) + DK;
        
        DWN.k2      = k2(ik1,imode);
        
        %approx de la solution par extrapolation & premiere approx de l intervale de recherche
        if ik>=4
            nextw0  = interp1(k2(ik-3:ik,imode),w0(ik-3:ik,imode),k2(ik1,imode),'pchip','extrap');
            nextwl  = interp1(k2(ik-3:ik,imode),w0(ik-3:ik,imode),k2(ik1,imode),'linear','extrap');
            dw1     = abs(nextw0-nextwl);
            dw0     = max(mean(err1(ik-3:ik,imode))*facDw*DK,dw1);
            dw0(dw0<1e-6)=1e-6;%pb precision machine
            
            wi      = nextw0-dw0;
            wf      = nextw0+dw0;
        elseif ik==3
            w1      = w0(1:3,imode);
            k21     = k2(1:3,imode);
            dwdk    = ((w1(2)-w1(3))*(k21(1)-k21(3))^2-(w1(1)-w1(3))*(k21(2)-k21(3))^2)/...
                ((k21(2)-k21(3))*(k21(1)-k21(3))^2-(k21(1)-k21(3))*(k21(2)-k21(3))^2);
            d2wdk2  = ((w1(2)-w1(3))-dwdk*(k21(2)-k21(3)))*2/(k21(2)-k21(3))^2;
            nextw0  = w1(3)+dwdk*(k2(ik1,imode)-k2(ik,imode))+d2wdk2*(k2(ik1,imode)-k2(ik,imode))^2/2;
            
            dw0     = err1(ik,imode)*facDw*DK;
            dw0(dw0<1e-6)=1e-6;%pb precision machine
            
            wi      = nextw0-dw0;
            wf      = min(nextw0+dw0,k2(ik1,imode)*vmax);
        elseif ik==2
            dwdk    = (w0(1,imode)-w0(2,imode))/(k2(1,imode)-k2(2,imode));
            nextw0  = w0(2,imode)+dwdk*(k2(ik1,imode)-k2(ik,imode));
            
            dw0     = abs(w0(ik,imode)-nextw0)/2;
            dwmin   = k2(ik1,imode)*(vmax-vmax2)/2;
            dw0     = min(dw0,dwmin);
            dw0(dw0<1e-6)=1e-6;%pb precision machine
            
            wi      = nextw0-dw0;
            wf      = min(nextw0+dw0,k2(ik1,imode)*vmax);
        elseif ik==1
            nextw0  = k2(ik1,imode)*vmax;
            
            dw0     = abs(w0(ik,imode)-nextw0)/2;
            dwmin   = k2(ik1,imode)*(vmax-vmax2)/2;
            dw0     = min(dw0,dwmin);
            dw0(dw0<1e-6)=1e-6;%pb precision machine
            
            wi  = nextw0-dw0;
            wf  = nextw0;
        end
        
        %         plotthepb_H(wi,wf,para,DWN)
        %verifier que la partie reelle et imaginaire change de signe
        [signp,sign0,sign1]=checkchgmtsign_H_vec(wi,wf,para,DWN,2);
        %         plot(k2(ik1,imode),nextw0,'c.')
        
        %         dw00=dw0;
        
        %verifier qu on n est pas proche d une branche asymptotique
        %ds ce cas la il y a une des parties qui croise 0 mais rechange
        %vite de signe, l intervalle de recherche doit etre tres petit
        %qd le sign==0 alors c est que c est des ondes homogenes
        if ((sign0>0 && sign1<=0 || sign0<=0 && sign1>0))
            %   plotthepb_H(wi,wf,para,DWN)               plotthepb_H(wi-dw0,wf+dw0,para,DWN)
            
            [~,~,~,det1,para1,DWN1]=checkchgmtsign_H_vec(wi,wf,para,DWN,50);
            solw    = cherche_zero_dic(para1,DWN1,det1);
            %   figure(203);hold on;plot(DWN1.omegac,real(det1),'k');plot(DWN1.omegac,imag(det1),'r');plot(DWN1.omegac,abs(det1),'')
            
            
            nind=length(solw);
            if nind~=0
                %on elimine les asymptotes de vP
                dwl=2e-5;
                if para.pol==2
                    while nind~=0
                        if (min(min(abs(solw*ones(1,ncapas)-ones(nind,1)*vP.'*k2(ik1,imode))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(solw*ones(1,ncapas)-ones(nind,1)*vP.'*k2(ik1,imode)),[],2);
                            [~,ind0]=min(abs(solw-vP(ind0)*k2(ik1,imode)));
                            solw(ind0)=[];
                            nind=length(solw);
                        else
                            break
                        end
                    end
                end
                %on elimine les asymptotes de vS
                if nind~=0
                    while nind~=0
                        if (min(min(abs(solw*ones(1,ncapas)-ones(nind,1)*vT.'*k2(ik1,imode))))<dwl)
                            %on est sur une asymptote, on enleve le point
                            [~,ind0]=min(abs(solw*ones(1,ncapas)-ones(nind,1)*vT.'*k2(ik1,imode)),[],2);
                            [~,ind0]=min(abs(solw-vT(ind0)*k2(ik1,imode)));
                            solw(ind0)=[];
                            nind=length(solw);
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
                    w0(ik1,imode)   = solw;
                    err             = abs(w0(ik1,imode)-nextw0)*2/DK;
                    err             = max(err,err1(ik,imode)/2);
                    err1(ik1,imode) = err;
                    ik              = ik+1;
                    continue
                end
                %                 facDw   = 1;
            end
        end
        
        %pb qd l intervale de recherche est trop petit
        j=1;
        %         dw0=dw00;
        %         while (sign0>0 || sign1>0) && j<4
        %             % plotthepb_H(wi,wf,para,DWN)
        %             dw0 = dw0*2;
        %             wi  = nextw0-dw0;
        %             wf  = nextw0+dw0;
        %             [signp,sign0,sign1,~,~,~]=checkchgmtsign_H_vec(wi,wf,para,DWN,2);
        %             j=j+1;
        %         end
        
        if signp<=0
            %tt est ok on recherche
            [~,~,~,det1,para1,DWN1] = checkchgmtsign_H_vec(wi,wf,para,DWN,3);
            solw                    = cherche_zero_dic(para1,DWN1,det1);
            nind                    = length(solw);
            if nind~=1
                DK              = DK/2;
                ik              = ik-1;
                err1(ik-1,imode)  = err1(ik,imode)/3;
                continue
            else
                w0(ik1,imode)   = solw;
                if j==1
                    err             = abs(w0(ik1,imode)-nextw0)*2/DK;
                    err             = max(err,err1(ik,imode)/2);
                    err1(ik1,imode) = err;
                    DK              = min(facav*DK,DK1);
                else
                    %on a trop reduit l intervale de recherche precedement
                    err             = dw0*2/DK;
                    err1(ik1,imode) = err;
                end
                ik              = ik+1;
            end
            %             facDw   = 1;
        else
            %si rien n a marche avant, on diminue le pas d avance
            DK      = DK/2;
%             if ik>5
%                 ik = ik-1;
%             end
            %             facDw   = 8;
        end
        % plot(k2(1:(ik1),imode),[w0(1:(ik),imode);nextw0],'r.')
        % plot(k2(1:ik,imode),w0(1:ik,imode),'r.')
        % plot(k2(ik1,imode),w0(ik1,imode),'r.')
    end
    % plot(k2(1:(ik),imode),w0(1:(ik),imode),'r')
    ikmax(imode)=ik1;
end
toc

figure(205); hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica asintotas                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2as  = linspace(0,kmax,1e3);
for ic=1:para.nsubmed
    C   = para.sub(ic).C;
    rho = para.sub(ic).rho;
    vT0  = sqrt(C(6,6)/rho);
    if  para.pol==1
        w=k2as*vT0;
        i0=find(w>wmax,1);
        if isempty(i0)
            i0=length(k2as);
        end
        plot(k2as(1:i0),w(1:i0),'k');
    else
        vL  = sqrt(C(1,1)/rho);
        w   = k2as*vT0;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2as);
        end
        plot(k2as(1:i0),w(1:i0),'k');
        w   = k2as*vL;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2as);
        end
        plot(k2as(1:i0),w(1:i0),'r')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           calculo de vg=dw/dk =f(w)                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vg=zeros(nmode,max(ikmax)-1);
f1=zeros(nmode,max(ikmax)-1);

for imode=1:nmode
    vg(imode,1:(ikmax(imode)-1))=diff(w0(1:ikmax(imode),imode))./diff(k2(1:ikmax(imode),imode));
    f1(imode,1:(ikmax(imode)-1))=(w0(1:ikmax(imode)-1,imode)+w0(2:ikmax(imode),imode))/2/2/pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(205);hold on;
for imode=1:nmode
    plot(k2(1:ikmax(imode),imode),w0(1:ikmax(imode),imode),'k')
end
plot(k20,w00,'.r')
vasin = unique(vT);
for iasin=2:length(vasin)-1
    [k20,w00]   = dispersion_curve_k_critik_Haskel(para,wmax,vasin(iasin));
    plot(k20,w00,'.r')
end

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
    plot( f1(imode,1:(ikmax(imode)-1)),vg(imode,1:(ikmax(imode)-1)),'k')
end