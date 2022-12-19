function [UXW,SXW]=calcul_US_DWN_PSV_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,tipo_onda)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)


ncapas      = para.nsubmed;
nk2         = length(kxk);
nrec        = length(xr);
nrecz       = length(zr0);

xs          = coordf.xs;
nxs       	= length(xs);
zs          = coordf.zs;

UXW         = zeros(2,nrec,nxs);
SXW         = zeros(2,2,nrec,nxs);

ksi         = para.reg(1).sub(ncapas).ksi;
kpi         = para.reg(1).sub(ncapas).kpi;

%*****************************************************%
% cas special des ondes de rayleigh ds un semi espace %
%*****************************************************%

if para.nsubmed==1 && tipo_onda==3
    %incidencia de ondas Rayleigh sin dispersion
    %expression analitica
    kri         = para.reg(1).kri;

    chi1    = sqrt(1-(kpi/kri)^2);
    chi2    = sqrt(1-(ksi/kri)^2);
    chi1    = sign(real(chi1)).*chi1;
    chi2    = sign(real(chi2)).*chi2;
    
    sq12    = sqrt(chi1*chi2);
    sq21    = sqrt(chi1/chi2);
    
    krix    = sign(kxk)*kri;
    kriz    = kri;
    
    u       = zeros(2,nk2);
    
    for irz=1:nrecz
        ixr     = find(izr0==irz);
        nrecx   = length(ixr);
        zricr   = zr0(irz);

        u(1,:)  =         (       exp(-kriz*chi1*abs(zricr))-sq12.*exp(-kriz*chi2*abs(zricr)))*sign(kxk);
        u(2,:)  = 1i*sq21*(-sq12.*exp(-kriz*chi1*abs(zricr))+      exp(-kriz*chi2*abs(zricr)));
        
        for irx=1:nrecx
            ir  = ixr(irx);
            expx=exp(-1i*(krix.*(xr(ir)-xs)));
            if salu(ir)==1
                UXW(1,ir,:)	= u(1,:).*expx;
                UXW(2,ir,:)	= u(2,:).*expx;
            end
        end
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul pour les autres incidences %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% identification des composantes du nombre d onde incident
if tipo_onda==1
    %incidencia de ondas P
    k2          = kxk*kpi;
elseif tipo_onda==2
    %incidencia de ondas S
    k2          = kxk*ksi;
elseif tipo_onda==3
    %incidencia de ondas Rayleigh, normali
    k2          = kxk*ksi;
end

%***************************************************************************************%
%  calcul des nombres d'onde verticaux, polarization & matrice de la solution homogene  %
%***************************************************************************************%

DWNinc.k2   = k2;
DWNinc      = calcul_A_DWN_PSV_Ncapas_HS(para,DWNinc);
A_DWN       = DWNinc.A_DWN;
k1          = DWNinc.k1 ;

%**********************************************************************************************%
% calcul des conditions a l'interface SE dues aux OP incidentes et resolution des coefficients %
%**********************************************************************************************%

MAT         = para.reg(1).sub;

% l incidence provient du semi-espace
ics         = ncapas;

Ci          = MAT(ics).Ci;
mic         =-(-1)^ics;%convention de signe

% position en z de la derniere interface sup au dessus du demi-espace
ztop= 0;
for inch=1:ics-1
    ztop=ztop+MAT(inch).h;
end

zt  = ztop-zs;% position relative dans la couche de la source

% calcul OP incidente, u et derivees de u pour calcul de t %
if tipo_onda==1
    %incidencia de ondas P
    k1f     = k1(2,:,ics);%k1P
    u1      =-DWNinc.u1(2,:,ics);
    u2      = DWNinc.u2(2,:,ics);%norma 1
elseif tipo_onda==2
    %incidencia de ondas S
    k1f     = k1(1,:,ics);%k1S
    u1      = DWNinc.u1(1,:,ics);
    u2      =-DWNinc.u2(1,:,ics);
elseif tipo_onda==3
    %### a completer
    %cf OPI
end
%propagation z
expt    = exp( 1i.*k1f.*zt);
%deplacement
u1      = u1.*expt;
u2      = u2.*expt;
%terme utile au calcul de s11 et s12
u11     = 1i*k1f.*u1;
u21     = 1i*k1f.*u2;
u12     =-1i*k2 .*u1;
u22     =-1i*k2 .*u2;

%indice de los terminos en el vector fuente
if ics==1
    is110   = 1;%SIGMA11 EN 0
    is120   = 2;%SIGMA12 EN 0
elseif ics>1
    is110   = (ics-2)*4+2+1;%SIGMA11 EN 0
    is120   = (ics-2)*4+2+2;%SIGMA12 EN 0
    iu10    = (ics-2)*4+2+3;%U1 EN 0
    iu20    = (ics-2)*4+2+4;%U2 EN 0
end

% terme source
Fsource         = zeros(4*(ncapas-1)+2,nk2);

Fsource(is110,:)= Ci(1,1)*u11 + Ci(1,2)*u22;                                 %sigma11 en 0
Fsource(is120,:)= Ci(6,6)*(u12+u21);                                         %sigma12 en 0
if ics>1
    Fsource(iu10,:)=u1;                                                     %U1 EN 0
    Fsource(iu20,:)=u2;                                                     %U2 EN 0
end
Fsource         = Fsource*mic;

%Resolution du systeme des conditions aux limites
% Solution   = A_DWN\Fsource;
Solution=zeros(4*(ncapas-1)+2,nk2);
for i=1:nk2
    Solution(:,i)   = A_DWN(:,:,i)\Fsource(:,i);
end

%******************************************%
% calcul des champs au niveau du recepteur %
%******************************************%

u1          = DWNinc.u1;
u2          = DWNinc.u2;
for irz=1:nrecz
    ixr     = find(izr0==irz);
    nrecx   = length(ixr);
    
    goU     = sum(salu(ixr(1:nrecx)));
    goS     = sum(sals(ixr(1:nrecx)));
    
    %Ecriture de la solution en fonction de la position du recepteur
    zricr   = zr0(irz);
    
    %determination de la couche du recepteur et de la profondeur relative
    icr     = 1;
    while zricr>MAT(icr).h && icr<ncapas
        zricr   = zricr-MAT(icr).h;
        icr     = icr+1;
    end
    
    exp1=exp(-1i.*k1(1,:,icr).* zricr);
    exp2=exp(-1i.*k1(2,:,icr).* zricr);
    if icr<ncapas
        exp3=exp( 1i.*k1(1,:,icr).*(zricr-MAT(icr).h));
        exp4=exp( 1i.*k1(2,:,icr).*(zricr-MAT(icr).h));
    end
    if goU>0
        %calcul des deplacements diffractes
        if icr<ncapas
            U1KW = ...
                +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ...
                +Solution(3+(icr-1)*4,:).*u1(1,:,icr).* exp3 ...
                -Solution(4+(icr-1)*4,:).*u1(2,:,icr).* exp4 ;
            U2KW = ...
                +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ...
                -Solution(3+(icr-1)*4,:).*u2(1,:,icr).* exp3 ...
                +Solution(4+(icr-1)*4,:).*u2(2,:,icr).* exp4;
        else
            U1KW = ...
                +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ;
            U2KW = ...
                +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ;
        end
        U1KW(isnan(U1KW))=0;
        U2KW(isnan(U2KW))=0;
    end
    
    if goS>0
        S11         = zeros(2,1);
        S12         = zeros(2,1);
        S22         = zeros(2,1);
        %calcul des contraintes diffractees
        for i=1:2
            S11(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,1).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,2).*u2(i,:,icr));
            S22(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,2).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,1).*u2(i,:,icr));
            S12(i,:)= -1i*MAT(icr).Ci(6,6)*( k1(i,:,icr).*u2(i,:,icr) +k2.*u1(i,:,icr));
        end
        if icr<ncapas
            S11KW = ...
                +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2 ...
                -Solution(3+(icr-1)*4,:).*S11(1,:).* exp3 ...
                +Solution(4+(icr-1)*4,:).*S11(2,:).* exp4;
            S22KW = ...
                +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2 ...
                -Solution(3+(icr-1)*4,:).*S22(1,:).* exp3 ...
                +Solution(4+(icr-1)*4,:).*S22(2,:).* exp4;
            S12KW = ...
                +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2 ...
                +Solution(3+(icr-1)*4,:).*S12(1,:).* exp3 ...
                -Solution(4+(icr-1)*4,:).*S12(2,:).* exp4;
        else
            S11KW = ...
                +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2;
            S22KW = ...
                +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2;
            S12KW = ...
                +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2;
        end
        S11KW(isnan(S11KW))=0;
        S22KW(isnan(S22KW))=0;
        S12KW(isnan(S12KW))=0;
    end
    
    %rajout des champs incidents
    if icr==ics
        zrs=zr0(irz)-zs;
         % calcul OP incidente, u et derivees de u pour calcul de t %
        if tipo_onda==1
            %incidencia de ondas P
            u1inc 	=-DWNinc.u1(2,:,ics);
            u2inc 	= DWNinc.u2(2,:,ics);%norma 1
        elseif tipo_onda==2
            %incidencia de ondas S
            u1inc  	= DWNinc.u1(1,:,ics);
            u2inc 	=-DWNinc.u2(1,:,ics);
        elseif tipo_onda==3
            %### a completer
            %cf OPI
        end
        %propagation z
        expt    = exp( 1i.*k1f.*zrs);
        %deplacement
        u1inc  	= u1inc.*expt;
        u2inc  	= u2inc.*expt;
        %terme utile au calcul de s11 et s12
        u11     = 1i*k1f.*u1inc;
        u21     = 1i*k1f.*u2inc;
        u12     =-1i*k2 .*u1inc;
        u22     =-1i*k2 .*u2inc;
        
        if goU>0
            U1KW    = U1KW  + u1inc;
            U2KW    = U2KW  + u2inc;
        end
        if goS>0
            S11KW   = S11KW + Ci(1,1)*u11 + Ci(1,2)*u22;
            S22KW   = S22KW + Ci(1,1)*u22 + Ci(1,2)*u11;
            S12KW   = S12KW + Ci(6,6)*(u12+u21);
        end
    end
    
    if min(xs==xs(1))==1 %todos son iguales
            for irx=1:nrecx
                ir = ixr(irx);
                expx=exp(-1i*(xr(ir)-xs).*k2);
                if salu(ir)==1
                    UXW(1,ir,:)	= U1KW.*expx;
                    UXW(2,ir,:)	= U2KW.*expx;
                end
                
                if sals(ir)==1
                    SXW(1,1,ir,:)	= S11KW.*expx;
                    SXW(2,2,ir,:)	= S22KW.*expx;
                    SXW(1,2,ir,:)	= S12KW.*expx;
                end
            end
    else
        for ixs=1:nxs
            for irx=1:nrecx
                ir = ixr(irx);
                expx=exp(-1i*(xr(ir)-xs(ixs))*k2(ixs));
                if salu(ir)==1
                    UXW(1,ir,ixs)	= U1KW.*expx;
                    UXW(2,ir,ixs)	= U2KW.*expx;
                end
                
                if sals(ir)==1
                    SXW(1,1,ir,ixs)	= S11KW.*expx;
                    SXW(2,2,ir,ixs)	= S22KW.*expx;
                    SXW(1,2,ir,ixs)	= S12KW.*expx;
                end
            end
        end
    end
end

%cuidado 1 y 2 son intercambiados entre IBEM y DWN
UXW(1:2,:,:)= UXW(2:-1:1,:,:);
SXW0        = SXW;
SXW(1,1,:,:)= SXW0(2,2,:,:);
SXW(2,2,:,:)= SXW0(1,1,:,:);