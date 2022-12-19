function dispersion_curve_kfix_Haskel(para)

% funcion que busqua a establecer las curvas de dispersion por un medio
% estratificado sobre un semi espacio
% por cada k(kfix),
% se busca los ceros de la fase de un sub-determinante de la matriz de Haskel 
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
tic
pol     = para.pol;
para.sub= para.reg(1).sub;
beta    = zeros(para.nsubmed,1);
for ms=1:para.nsubmed
    para.sub(ms).C(6,6)  = para.sub(ms).rho*para.sub(ms).bet^2;
    beta(ms)  = para.reg(1).sub(ms).bet;
end
vmin    = min(beta)*(1+1e-6);
vmax    = max(beta)*(1-1e-6);
vdeep   = beta(end)*(1-1e-6);

wmax        = 2*pi*para.fmax;

nk          = 200;
kmax        = wmax/min(beta);
DK          = kmax/nk;

[k2c,w0c]   = dispersion_curve_k_critik_Haskel(para,wmax,vdeep);
nmode       = length(k2c);
wsol      	= zeros(nk,nmode);
ikmax       = zeros(1,nmode);

nf          = 1e2*nmode;
indk        = 1:(nf-1);

for j=1:nk
    k20         = j*DK;
    DWN0        = struct('k2',k20);
    DWN0.omegac	= linspace(k20*vmin*.5,min(k20*vmax,wmax),nf);
    
    if  pol==1
        tmp=mode_Love(para,DWN0);
    else
        tmp=mode_Rayleigh_2(para,DWN0);
    end
    tmp(tmp==inf)=0;
    tmp1    = angle(tmp)-pi/2;
    indd1   = logical(tmp1(1:(nf-1)).*tmp1(2:nf)<=0);
    indd1   = indk(indd1);
    if ~isempty(indd1)
        w0=cherche_zero(DWN0.omegac,real(tmp).',indd1).';
        w0(w0==0)=[];
        nind    = length(w0);
        if nind>nmode
            w0=w0(1:nmode);
        end
        wsol(j,:)=[w0;zeros(nmode-nind,1)];
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2  = (1:nk)*DK;
figure(205);hold on;
for j=1:nmode
    indi    = find(wsol(:,j)~=0,1,'first');
    indf    = find(wsol(:,j)~=0,1,'last');
    ikmax(j)= indf-indi+1;

    plot([k2c(j),k2(indi:indf)],[w0c(j),wsol(indi:indf,j).'],'r')
end
xlabel('k_x');
ylabel('\omega');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica asintotas                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2  = linspace(0,kmax,1e3);
for ic=1:para.nsubmed
    C   = para.sub(ic).C;
    rho = para.sub(ic).rho;
    vT  = sqrt(C(6,6)/rho);
    if  para.pol==1
        w=k2*vT;
        i0=find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'k');
    else
        vL  = sqrt(C(1,1)/rho);
        w   = k2*vT;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'k');
        w   = k2*vL;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'r')
    end
end


k2  = (1:nk)*DK;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vp=w/k =f(w)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(206);hold on;
for imode=1:nmode
    indi=find(wsol(:,imode)~=0,1,'first');
    indf=find(wsol(:,imode)~=0,1,'last');
    vp=wsol(indi:indf,imode)./k2(indi:indf).';
    if imode==1
        vp(1)=vdeep;
    end
    plot(wsol(indi:indf,imode)/2/pi,vp,'k')
end
xlabel('frecuencia');
ylabel('Vp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vg=zeros(nmode,max(ikmax)-1);
f1=zeros(nmode,max(ikmax)-1);

for imode=1:nmode
    indi=find(wsol(:,imode)~=0,1,'first');
    indf=find(wsol(:,imode)~=0,1,'last');
    vg(imode,1:indf-indi)=diff(wsol(indi:indf,imode))./diff(k2(indi:indf)).';
    f1(imode,1:indf-indi)=(wsol(indi:(indf-1),imode)+wsol((indi+1):indf,imode))/2/2/pi;
end
figure(207);hold on;
for imode=1:nmode
    plot( f1(imode,1:(ikmax(imode)-1)),vg(imode,1:(ikmax(imode)-1)),'k')
end
xlabel('frecuencia');
ylabel('Vg');
