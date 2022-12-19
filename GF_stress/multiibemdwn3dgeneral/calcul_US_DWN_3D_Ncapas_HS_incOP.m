function [UXW,SXW]=calcul_US_DWN_3D_Ncapas_HS_incOP(para,xr0,yr0,zr0,izr0,salu,sals,coordf,kxk,tipo_onda,kzsigno,phi)
% Campo en medio estratificado por incidencia de onda de presión o de corte
% UXW (3,n) desplazamiento
% SXW (3,3,n) esfuerzos

% Procedimiento:
% Se 'rotan' las posiciones de los receptores, 
% Se calcula en el plano de polarización de la onda plana
% Se rota el resultado en los receptores al ángulo inicial.
%
% xr,yr,zr0  coordenadas de los receptores a la profundidad zr0
% izr0  índice de profundidad de los receptores
% salu, sals  indicadores de variables de salida solicitadas
% coordf  coordenadas de la fuente actual
% kxk  cartesiana x de la normal de la OP
% kyk  cartesiana y de la normal de la OP
% tipoonda  polarización de la onda plana

% rotar coordenadas
if phi ~= 0
[xr,~] = rotarXY(xr0,yr0,-phi);
else
xr = xr0;
end
ncapas      = para.nsubmed;
nk2         = length(kxk);
nrec        = length(xr);
nrecz       = length(zr0);

xs          = coordf.xs;
ys          = coordf.ys;
zs          = coordf.zs;
if (xs ~= 0 || ys ~= 0)
    xs = 0; %ys = 0;
end
nxs       	= length(xs);

UXW         = zeros(3,nrec,nxs);
SXW         = zeros(3,3,nrec,nxs);

ksi         = para.reg(1).sub(ncapas).ksi; % omega/beta
kpi         = para.reg(1).sub(ncapas).kpi; % omega/alpha

%% *****************************************************%
% caso especial de onda de Rayleigh en un semiespacio %
%*****************************************************%
if para.nsubmed==1 && tipo_onda==4
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
    kr          = kxk*kpi;
elseif tipo_onda==2
    %incidencia de ondas SV
    kr          = kxk*ksi;
elseif tipo_onda==3
    %incidencia de ondas SH
    kr          = kxk*ksi;
else
    %incidencia de ondas Rayleigh, normali
    kr          = kxk*ksi;
end

%***************************************************************************************%
%  calcul des nombres d'onde verticaux, polarization & matrice de la solution homogene  %
%***************************************************************************************%

DWNinc.k2   = kr;
% DWN.omegac  = 2*pi*fj;

% Entra  DWN.kr, para.nsubmed, para.reg
% Sale   DWN.kz, DWN.u1, DWN.u2, DWN.A_DWN (P-SV), DWN.B_DWN (SH)
%DWNinc      = calcul_A_DWN_3D_polar_Ncapas_HS(para,DWNinc);
if tipo_onda~=3 % [inplane]
DWNinc      = calcul_A_DWN_PSV_Ncapas_HS(para,DWNinc);
A_DWN       = DWNinc.A_DWN; % matriz P_SV
kz          = DWNinc.k1 ;   % componente vertical
else
B_DWN       = DWNinc.B_DWN; % matriz SH
kz          = DWNinc.k1 ;   % componente vertical
end
%**********************************************************************************************%
% calcul des conditions a l'interface SE dues aux OP incidentes et resolution des coefficients %
%**********************************************************************************************%

% l incidence provient du semi-espace
MAT         = para.reg(1).sub; %material del fondo
ics         = ncapas; % cantidad de capas del fondo estratificado
Ci          = MAT(ics).Ci; %(6,6) : mu
mic         =-(-1)^ics;%convention de signe

% position en z de la derniere interface sup au dessus du demi-espace
ztop= 0; %profundidad de la interfaz del semiespacio
for inch=1:ics-1
    ztop=ztop+MAT(inch).h;
end
zt  = ztop-zs;% position relative dans la couche de la source

% calcul OP incidente, u et derivees de u pour calcul de t %
if tipo_onda==1
    %incidencia de ondas P
    k1f     = kz(2,:,ics);%k1P
    ur      =-DWNinc.u1(2,:,ics); % desp en r
    uz      = DWNinc.u2(2,:,ics); % desp en z   norma 1
elseif tipo_onda==2
    %incidencia de ondas SV
    k1f     = kz(1,:,ics);%k1S
    ur      = DWNinc.u1(1,:,ics); % desp en r
    uz      =-DWNinc.u2(1,:,ics); % desp en z
elseif tipo_onda==3
    %incidencia de ondas SH
    k1f	= kzsigno.*sqrt(ksi^2-kr.^2); %kz = + o -, la 2 soluciones son correctas
    exp2t   = exp(1i.*k1f.*zt);
    G220    = exp2t;%pair en k2
    % derivee p/x1 en X (x1=z)
    %pair
    ut=1i*k1f.*exp2t;
    
    k1 = zeros(ncapas,nk2);
    for ic=1:ncapas
        ksi     = para.reg(1).sub(ic).ksi;
        k1(ic,:)  = sqrt(ksi^2-kr.^2); %k1S
    end
    k1 =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));
else
    %### a completer
    %cf OPI
    error('falta programar la incidencia de ondas de Rayleigh')
end
%propagation z
expt    = exp( 1i.*k1f.*zt);

if tipo_onda~=3 % [inplane]
    %deplacement [con fase vertical]
    ur      = ur.*expt;
    uz      = uz.*expt;
    %terme utile au calcul de s11 et s12 [derivees]
    u11     = 1i*k1f.*ur;
    u21     = 1i*k1f.*uz;
    u12     =-1i*kr .*ur;
    u22     =-1i*kr .*uz;
    
    %indices de los terminos en el vector fuente
    % [La incidencia se coloca en la interfaz del semiespacio]
    if ics==1 % no hay capas (sólo tracciones)
        is110   = 1;%SIGMA11 EN 0
        is120   = 2;%SIGMA12 EN 0
        is130   = 1;
    elseif ics>1 % hay un estrato arriba
        is110   = (ics-2)*4+2+1;%SIGMA11 EN 0
        is120   = (ics-2)*4+2+2;%SIGMA12 EN 0
        is130   = (ics-2)*2+1+1;
        iu10    = (ics-2)*4+2+3;%U1 EN 0
        iu20    = (ics-2)*4+2+4;%U2 EN 0
        iu30    = (ics-2)*2+1+2;
    end
    
    FsourcePSV  = zeros(4*(ncapas-1)+2,nk2);
    SolutionPSV = zeros(4*(ncapas-1)+2,nk2);
    
    % terme source P-SV
    % (tracciones)
    FsourcePSV(is110,:)= Ci(1,1)*u11 + Ci(1,2)*u22;                                 %sigma11 en 0
    FsourcePSV(is120,:)= Ci(6,6)*(u12+u21);                                         %sigma12 en 0
    if ics>1 % (desplazamientos)
        FsourcePSV(iu10,:)=ur;                                                     %U1 EN 0
        FsourcePSV(iu20,:)=uz;                                                     %U2 EN 0
    end
    FsourcePSV         = FsourcePSV*mic;
    %Resolution du systeme des conditions aux limites
    % Solution   = A_DWN\Fsource;
    for i=1:nk2
        SolutionPSV(:,i)   = A_DWN(:,:,i)\FsourcePSV(:,i);
    end
    
    ur          = DWNinc.u1; % Campo incidente de desplazamiento
    uz          = DWNinc.u2; % en la interfaz y sin fase vertical.
end

if tipo_onda==3 % [anti-plane]
    FsourceSH   = zeros(2*(ncapas-1)+1,nk2);
    SolutionSH  = zeros(2*(ncapas-1)+1,nk2);
    % termino de fuente SH
    FsourceSH(is130,:)    = Ci(6,6)*ut;   %sigmazr en 0 %falsos indices
    if ics>1
        FsourceSH(iu30,:) = G220;   	%Ut en 0 %falsos indices
    end
    FsourceSH= FsourceSH*mic;
    for i=1:nk2
        SolutionSH(:,i)   = B_DWN(:,:,i)\FsourceSH(:,i);
    end
end
%******************************************%
% calcul des champs au niveau du recepteur %
%******************************************%


for irz=1:nrecz % para cada profundidad de receptores
    ixr     = find(izr0==irz); %indices de los receptores a esta prof.
    nrecx   = length(ixr);
    
    goU     = sum(salu(ixr(1:nrecx))); % > 0 si se calcula U
    goS     = sum(sals(ixr(1:nrecx))); %                   S
    
    %Ecriture de la solution en fonction de la position du recepteur
    zricr   = zr0(irz);
    
    %determination de la couche du recepteur et de la profondeur relative
    icr     = 1;
    while zricr>MAT(icr).h && icr<ncapas
        zricr   = zricr-MAT(icr).h;
        icr     = icr+1; %[couche]
    end
    
    % propagación vertical
    exp1=exp(-1i.*kz(1,:,icr).* zricr); %S desde arriba
    exp2=exp(-1i.*kz(2,:,icr).* zricr); %P desde arriba
    if icr<ncapas
        exp3=exp( 1i.*kz(1,:,icr).*(zricr-MAT(icr).h)); %S desde abajo
        exp4=exp( 1i.*kz(2,:,icr).*(zricr-MAT(icr).h)); %P desde abajo
    end
    
    if goU>0
        %calcul des deplacements diffractes
        if tipo_onda~=3 % [inplane]
            if icr<ncapas
                U1KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*ur(1,:,icr).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*ur(2,:,icr).* exp2 ...
                    +SolutionPSV(3+(icr-1)*4,:).*ur(1,:,icr).* exp3 ...
                    -SolutionPSV(4+(icr-1)*4,:).*ur(2,:,icr).* exp4 ;
                U2KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*uz(1,:,icr).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*uz(2,:,icr).* exp2 ...
                    -SolutionPSV(3+(icr-1)*4,:).*uz(1,:,icr).* exp3 ...
                    +SolutionPSV(4+(icr-1)*4,:).*uz(2,:,icr).* exp4;
            else
                U1KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*ur(1,:,icr).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*ur(2,:,icr).* exp2 ;
                U2KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*uz(1,:,icr).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*uz(2,:,icr).* exp2 ;
            end
            U1KW(isnan(U1KW))=0;
            U2KW(isnan(U2KW))=0;
        else % [anti-plane]
            if icr<ncapas
                UyKW = ...
                    +SolutionSH(1+(icr-1)*2,:).* exp1.' ...
                    +SolutionSH(2+(icr-1)*2,:).* exp3.';
            else
                UyKW = SolutionSH(1+(icr-1)*2,:).* exp1.';
            end
            UyKW(isnan(UyKW))=0;
        end
    end
    
    if goS>0
        if tipo_onda~=3 % [inplane]
            S11         = zeros(2,1);
            S12         = zeros(2,1);
            S22         = zeros(2,1);
            %calcul des contraintes diffractees
            for i=1:2
                S11(i,:)= -1i*(kz(i,:,icr).*MAT(icr).Ci(1,1).*ur(i,:,icr) +kr.*MAT(icr).Ci(1,2).*uz(i,:,icr));
                S22(i,:)= -1i*(kz(i,:,icr).*MAT(icr).Ci(1,2).*ur(i,:,icr) +kr.*MAT(icr).Ci(1,1).*uz(i,:,icr));
                S12(i,:)= -1i*MAT(icr).Ci(6,6)*( kz(i,:,icr).*uz(i,:,icr) +kr.*ur(i,:,icr));
            end
            if icr<ncapas
                S11KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S11(2,:).* exp2 ...
                    -SolutionPSV(3+(icr-1)*4,:).*S11(1,:).* exp3 ...
                    +SolutionPSV(4+(icr-1)*4,:).*S11(2,:).* exp4;
                S22KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S22(2,:).* exp2 ...
                    -SolutionPSV(3+(icr-1)*4,:).*S22(1,:).* exp3 ...
                    +SolutionPSV(4+(icr-1)*4,:).*S22(2,:).* exp4;
                S12KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S12(2,:).* exp2 ...
                    +SolutionPSV(3+(icr-1)*4,:).*S12(1,:).* exp3 ...
                    -SolutionPSV(4+(icr-1)*4,:).*S12(2,:).* exp4;
            else
                S11KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S11(2,:).* exp2;
                S22KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S22(2,:).* exp2;
                S12KW = ...
                    +SolutionPSV(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +SolutionPSV(2+(icr-1)*4,:).*S12(2,:).* exp2;
            end
            S11KW(isnan(S11KW))=0;
            S22KW(isnan(S22KW))=0;
            S12KW(isnan(S12KW))=0;
        else % [anti-plane]
            C0= -1i*MAT(icr).Ci(6,6)*k1(icr,:).';
            if icr<ncapas
                SzyKW = ...
                    +SolutionSH(1+(icr-1)*2,:).*C0.* exp1.' ...
                    -SolutionSH(2+(icr-1)*2,:).*C0.* exp3.';
            else
                SzyKW =SolutionSH(1+(icr-1)*2).*C0.* exp1.';
            end
            SzyKW(isnan(SzyKW))=0;
            
            %impair en k2
            C0= -1i*MAT(icr).Ci(6,6)*k2;
            if icr<ncapas
                SxyKW = ...
                    +SolutionSH(1+(icr-1)*2,:).*C0.* exp1.' ...
                    +SolutionSH(2+(icr-1)*2,:).*C0.* exp3.';
            else
                SxyKW =SolutionSH(1+(icr-1)*2).*C0.* exp1.';
            end
            SxyKW(isnan(SxyKW))=0;
        end
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
            %incidencia de ondas SV
            u1inc  	= DWNinc.u1(1,:,ics);
            u2inc 	=-DWNinc.u2(1,:,ics);
        elseif tipo_onda==3
            % SH
            u3inc   = 1;
        else
            %### a completer
            %cf OPI
        end
        %propagation z
        expt    = exp( 1i.*k1f.*zrs);
        
        if tipo_onda~=3 % [inplane]
            %deplacement
            u1inc  	= u1inc.*expt;
            u2inc  	= u2inc.*expt;
            %terme utile au calcul de s11 et s12
            u11     = 1i*k1f.*u1inc;
            u21     = 1i*k1f.*u2inc;
            u12     =-1i*kr .*u1inc;
            u22     =-1i*kr .*u2inc;
            
            if goU>0
                U1KW    = U1KW  + u1inc;
                U2KW    = U2KW  + u2inc;
            end
            if goS>0
                S11KW   = S11KW + Ci(1,1)*u11 + Ci(1,2)*u22;
                S22KW   = S22KW + Ci(1,1)*u22 + Ci(1,2)*u11;
                S12KW   = S12KW + Ci(6,6)*(u12+u21);
            end
        else % [anti plane]
            if goU>0
                G22r = u3inc.*expt;
                UyKW = UyKW + G22r;
            end
            
            if goS>0
                %%%derivee p/x1 en X (x1=z)
                %pair
                G22z=1i*k1f.*expt.';
                
                % derivee p/x2 en X (x2=x)
                %impair
                G22x=-1i*kr.*expt.';
                
                %calcul des tractions t1j
                SzyKW = SzyKW + MAT(icr).Ci(6,6)*G22z;
                SxyKW = SxyKW + MAT(icr).Ci(6,6)*G22x;
            end
        end
    end
    
    if min(xs==xs(1))==1 %todos son iguales
        for irx=1:nrecx
            ir = ixr(irx);
            expx=exp(-1i*(xr(ir)-xs).*kr);
            if salu(ir)==1
                if tipo_onda~=3
                UXW(3,ir,:)	= U1KW.*expx;
                UXW(1,ir,:)	= U2KW.*expx;
                else
                UXW(2,ir,:) = UyKW.*expx;
                end
            end
            
            if sals(ir)==1
                if tipo_onda~=3
                SXW(3,3,ir,:)	= S11KW.*expx;
                SXW(1,1,ir,:)	= S22KW.*expx;
                SXW(1,3,ir,:)	= S12KW.*expx;
                 SXW(3,1,ir,:)	= S12KW.*expx;
                else
                SXW(1,2,ir,:)   = SxyKW.*expx;
                 SXW(2,1,ir,:)  = SxyKW.*expx;
                SXW(3,2,ir,:)   = SzyKW.*expx;
                 SXW(2,3,ir,:)  = SzyKW.*expx;
                end
            end
        end
    else
        for ixs=1:nxs
            for irx=1:nrecx
                ir = ixr(irx);
                expx=exp(-1i*(xr(ir)-xs(ixs))*kr(ixs));
                if salu(ir)==1
                    if tipo_onda~=3
                    UXW(3,ir,ixs)	= U1KW.*expx;
                    UXW(1,ir,ixs)	= U2KW.*expx;
                    else
                    UXW(2,ir,ixs)   = UyKW.*expx;    
                    end
                end
                
                if sals(ir)==1
                    if tipo_onda~=3
                    SXW(3,3,ir,ixs)	= S11KW.*expx;
                    SXW(1,1,ir,ixs)	= S22KW.*expx;
                    SXW(1,3,ir,ixs)	= S12KW.*expx;
                     SXW(3,1,ir,ixs)= S12KW.*expx;
                    else
                    SXW(1,2,ir,ixs) = SxyKW.*expx;
                     SXW(2,1,ir,ixs)= SxyKW.*expx;
                    SXW(3,2,ir,ixs) = SzyKW.*expx;
                     SXW(2,3,ir,ixs)= SzyKW.*expx;
                    end
                end
            end
        end
    end
end

% En esta versión ya se ha coregido el siguiente:
% % cuidado 1 y 2 son intercambiados entre IBEM y DWN
% % UXW(1:2,:,:)= UXW(2:-1:1,:,:);
% % SXW0        = SXW;
% % SXW(1,1,:,:)= SXW0(2,2,:,:);
% % SXW(2,2,:,:)= SXW0(1,1,:,:);
% % Y los canales correctos para el 3D:
% % 2 es y
% % 3 es z

% rotar de regreso
if phi ~= 0
   if min(xs==xs(1))==1 %todos son iguales
        for irx=1:nrecx
            ir = ixr(irx);
            if salu(ir)==1
                [UXW(1,ir,:),UXW(2,ir,:)] = ...
                    rotarXY(UXW(1,ir,:),UXW(2,ir,:),phi);
            end
            
            if sals(ir)==1
                for ixs = 1:nxs
                    SXW(:,:,ir,ixs) = rotarS(SXW(:,:,ir,ixs),phi);
                end
            end
        end
    else
        for ixs=1:nxs
            for irx=1:nrecx
                ir = ixr(irx);
                if salu(ir)==1
                [UXW(1,ir,ixs),UXW(2,ir,ixs)] = ...
                    rotarXY(UXW(1,ir,ixs),UXW(2,ir,ixs),phi);
                end
                
                if sals(ir)==1
                    SXW(:,:,ir,ixs) = rotarS(SXW(:,:,ir,ixs),phi);
                end
            end
        end
    end
end
end

function [xp,yp] = rotarXY(x,y,phi)
phi = phi * pi / 180;
b11 = cos(phi);  b12 = sin(phi);
b21 = -b12;  b22 = b11;

xp = b11.* x + b12.* y;
yp = b21.* x + b22.* y;
end

function [Sp] = rotarS(S,phi)
phi = phi * pi / 180;
a = [[  cos(phi) sin(phi) 0]; ...
     [ -sin(phi) cos(phi) 0]; ...
     [   0        0       1]];
Sp = a * S * a.';
end
