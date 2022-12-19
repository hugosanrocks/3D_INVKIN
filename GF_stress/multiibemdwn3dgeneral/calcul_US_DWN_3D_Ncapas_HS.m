function [UXW1,SXW1,UXW2,SXW2,UXW3,SXW3]=calcul_US_DWN_3D_Ncapas_HS(para,xr,yr,zr0,izr0,salu,sals,coordf,fij,DWN)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% attention, pour calculer les differentes composantes
% il est nécessaire de subdiviser le calcul en fonction
% des parties paire et impaire (k2)
% cf G11, G22 paires et G12 impaire
% ce qui peut être fait en traitant séparement fint1 et fint2
ncapas      = para.nsubmed;
k2          = DWN.k2;
k1          = DWN.k1;
u1          = DWN.u1;
u2          = DWN.u2;

A_DWN       = DWN.A_DWN;
nk2       	= length(k2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion de la matrice %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nk2
    A_DWN(:,:,i)=inv(A_DWN(:,:,i));
end

nk2t        = 2*(nk2-1);%para.DWNnbptkx;
MAT         = para.reg(1).sub;
nrec        = length(xr);
nrecz       = length(zr0);

omegac      = DWN.omegac;

Kpos2       = k2;
Kneg2       = k2(2:(nk2-1));
DK          = DWN.dk2;
DK          = [DK,DK(end-1:-1:2)];

U1KW        = zeros(1,nk2);
U2KW        = zeros(1,nk2);
S11KW       = zeros(1,nk2);
S22KW       = zeros(1,nk2);
S12KW       = zeros(1,nk2);

U1KW0       = zeros(1,nk2t);
U2KW0       = zeros(1,nk2t);
S11KW0      = zeros(1,nk2t);
S22KW0      = zeros(1,nk2t);
S12KW0      = zeros(1,nk2t);

xs          = coordf.xs;
nxs       	= length(xs);
zs          = coordf.zs;

UXW1        = zeros(2,nrec,nxs);
SXW1        = zeros(2,2,nrec,nxs);
UXW2        = zeros(2,nrec,nxs);
SXW2        = zeros(2,2,nrec,nxs);

S11         = zeros(2,nk2);
S12         = zeros(2,nk2);
S22         = zeros(2,nk2);
Fsource     = zeros(4*(ncapas-1)+2,nk2);
Solution    = zeros(4*(ncapas-1)+2,nk2);
if abs(fij(2))<1e-3*abs(fij(1))
    fij(2)=0;
elseif abs(fij(1))<1e-3*abs(fij(2))
    fij(1)=0;
end
fijIBEM     = fliplr(fij);

%identification de la couche dans laquelle se trouve la force
ics =1;
hc  =MAT(1).h;
while zs>hc && ics<ncapas
    ics=ics+1;
    hc=hc+MAT(ics).h;
end
k1f(1,:)= k1(1,:,ics);%k1S
k1f(2,:)= k1(2,:,ics);%k1P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul du vecteur source et resolution des coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% le calcul du vecteur source se fait en utilisant
% la decomposition en ondes planes
if ics>1
    htop=0;
    for inch=1:ics-1
        htop=htop+MAT(inch).h;
    end
else
    htop=0;
end
hbot=0;
for inch=1:ics
    hbot=hbot+MAT(inch).h;
end

mxspecs1=1;
mxspecb1=1;
mxspecs2=1;
mxspecb2=1;

if nxs==1
    mxspec=1;
    if isfield(coordf,'vnx')
        lsegx=coordf.dr*coordf.vnz;
        lsegz=coordf.dr*coordf.vnx;
        mxspecs1=sinc((-lsegx*Kpos2-lsegz*k1f(1,:))/2/pi);
        mxspecb1=sinc((-lsegx*Kpos2+lsegz*k1f(1,:))/2/pi);
        mxspecs2=sinc((-lsegx*Kpos2-lsegz*k1f(2,:))/2/pi);
        mxspecb2=sinc((-lsegx*Kpos2+lsegz*k1f(2,:))/2/pi);
    end
end

%%%%%%%%%%%%%%%%%
% interface sup %
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des deplacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt      = htop-zs;%x=z dans l epaisseur
exp1t   = exp(-1i.*k1f(1,:).*abs(zt)).*mxspecs1;
exp2t   = exp(-1i.*k1f(2,:).*abs(zt)).*mxspecs2;
fac     = 1/(4*pi*MAT(ics).rho*omegac^2);
%pair en k2
G220=-1i*fac*( k1f(1,:)         .*exp1t	+ k2.^2./k1f(2,:)   .*exp2t);
%pair en k2
G110=-1i*fac*( k2.^2./k1f(1,:)	.*exp1t + k1f(2,:)          .*exp2t);
%impair en k2
G120=-1i*fac*k2*sign(zt).*(exp2t-exp1t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des derivees utiles au calcul des contraintes normales %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivee p/x1 en X (x1=z)
%pair
G2210=-sign(zt)*fac*(k1f(1,:).^2.*exp1t + k2.^2             .*exp2t);
%pair
G1110=-sign(zt)*fac*(k2.^2      .*exp1t + k1f(2,:).^2       .*exp2t);
%impair
G1210=-fac*k2.*(k1f(2,:).*exp2t-k1f(1,:).*exp1t);

% derivee p/x2 en X (x2=x)
%impair
G2220=-k2.*fac.*( k1f(1,:)      .*exp1t + k2.^2./k1f(2,:)   .*exp2t);
%impair
G1120=-k2.*fac.*(k2.^2./k1f(1,:).*exp1t + k1f(2,:)          .*exp2t);
%pair
G1220=-k2.^2*sign(zt).*fac.*(exp2t-exp1t);

%%%%%%%%%%%%%%%%%
% interface inf %
%%%%%%%%%%%%%%%%%
if ics~=ncapas
    zb      = hbot-zs;%x=z dans l epaisseur
    exp1b   = exp(-1i.*k1f(1,:).*abs(zb)).*mxspecb1;
    exp2b   = exp(-1i.*k1f(2,:).*abs(zb)).*mxspecb2;
    
    %pair en k2
    G22h=-1i*fac*( k1f(1,:)         .*exp1b + k2.^2./k1f(2,:)  	.*exp2b);
    %pair en k2
    G11h=-1i*fac*( k2.^2./k1f(1,:)  .*exp1b + k1f(2,:)          .*exp2b);
    %impair en k2
    G12h=-1i*fac*k2*sign(zb).*(exp2b-exp1b);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des derivees utiles au calcul des contraintes normales %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivee p/x1 en X (x1=z)
    %pair
    G221h=-sign(zb)*fac*(k1f(1,:).^2.*exp1b + k2.^2             .*exp2b);
    %pair
    G111h=-sign(zb)*fac*(k2.^2      .*exp1b + k1f(2,:).^2       .*exp2b);
    %impair
    G121h=-fac*k2.*(k1f(2,:).*exp2b-k1f(1,:).*exp1b);
    
    % derivee p/x2 en X (x2=x)
    %impair
    G222h=-k2.*fac.*( k1f(1,:)      .*exp1b + k2.^2./k1f(2,:)   .*exp2b);
    %impair
    G112h=-k2.*fac.*(k2.^2./k1f(1,:).*exp1b + k1f(2,:)          .*exp2b);
    %pair
    G122h=-k2.^2*sign(zb).*fac.*(exp2b-exp1b);
end

%indice de los terminos en el vector fuente
if ics==1
    is110=1;%SIGMA11 EN 0
    is120=2;%SIGMA12 EN 0
    if ncapas>1
        is11h=3;%SIGMA11 EN h(ics)
        is12h=4;%SIGMA12 EN h(ics)
        iu1h =5;%U1 EN h(ics)
        iu2h =6;%U2 EN h(ics)
    end
elseif ics>1
    is110=(ics-2)*4+2+1;%SIGMA11 EN 0
    is120=(ics-2)*4+2+2;%SIGMA12 EN 0
    iu10=(ics-2)*4+2+3;%U1 EN 0
    iu20=(ics-2)*4+2+4;%U2 EN 0
    if ics<=(ncapas-1)
        is11h=(ics-2)*4+2+5;%SIGMA11 EN h(ics)
        is12h=(ics-2)*4+2+6;%SIGMA12 EN h(ics)
        iu1h =(ics-2)*4+2+7;%U1 EN h(ics)
        iu2h =(ics-2)*4+2+8;%U2 EN h(ics)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci      = MAT(ics).Ci;
mic     =-(-1)^ics;


%fint1=fij(1);
%fint2=0;
if fij(1)~=0
    Fsource(is110,:)=Ci(1,1)*G1110 + Ci(1,2)*G1220;      %sigma11 en 0
    Fsource(is120,:)=Ci(6,6)*(G1120+G1210);                            %sigma12 en 0
    if ics>1
        Fsource(iu10,:)=G110;%U1 EN 0
        Fsource(iu20,:)=G120;%U2 EN 0
    end
    if (ncapas>1 && ics==1) || (ics<=(ncapas-1))
        Fsource(is11h,:)=Ci(1,1)*G111h + Ci(1,2)*G122h;	%sigma11 en h
        Fsource(is12h,:)=Ci(6,6)*(G112h+G121h);                        %sigma12 en h
        Fsource(iu1h,:) =G11h;%U1 EN h
        Fsource(iu2h,:) =G12h;%U2 EN h
    end
    Fsource=fij(1)*Fsource*mic;
    
    %Resolution du systeme des conditions aux limites
    %         for i=1:nk2
    %             Solution(:,i)   = squeeze(A_DWN(:,:,i))*squeeze(Fsource(:,i));
    %         end
    for i=1:4*(ncapas-1)+2
        Solution(i,:)   = sum(squeeze(A_DWN(i,:,:)).*Fsource,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des champs au niveau du recepteur %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for irz=1:nrecz
        ixr     = find(izr0==irz);
        nrecx   = length(ixr);
        
        goU=sum(salu(ixr(1:nrecx)));
        goS=sum(sals(ixr(1:nrecx)));
        %Ecriture de la solution en fonction de la position du recepteur
        zricr   = zr0(irz);
        
        icr     = 1;
        while zricr>MAT(icr).h && icr<ncapas
            zricr   = zricr-MAT(icr).h;
            icr     = icr+1;
        end
        
        exp1=exp(-1i.*k1(1,:,icr).* zricr);
        exp2=exp(-1i.*k1(2,:,icr).* zricr);
        exp3=exp( 1i.*k1(1,:,icr).*(zricr-MAT(icr).h));
        exp4=exp( 1i.*k1(2,:,icr).*(zricr-MAT(icr).h));
        
        if goU>0
            %calcul des deplacements diffractes
            if icr<ncapas
                U1KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ...
                    +Solution(3+(icr-1)*4,:).*u1(1,:,icr).* exp3 ...
                    -Solution(4+(icr-1)*4,:).*u1(2,:,icr).* exp4 ;
                U2KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*u2(1,:,icr).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*u2(2,:,icr).* exp4;
            else
                U1KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ;
                U2KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ;
            end
            U1KW(isnan(U1KW))=0;
            U2KW(isnan(U2KW))=0;
        end
        
        if goS>0
            %calcul des contraintes diffractees
            for i=1:2
                S11(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,1).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,2).*u2(i,:,icr));
                S22(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,2).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,1).*u2(i,:,icr));
                S12(i,:)= -1i*MAT(icr).Ci(6,6)*( k1(i,:,icr).*u2(i,:,icr) +k2.*u1(i,:,icr));
            end
            if icr<ncapas
                S11KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*S11(1,:).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*S11(2,:).* exp4;
                S22KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*S22(1,:).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*S22(2,:).* exp4;
                S12KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2 ...
                    +Solution(3+(icr-1)*4,:).*S12(1,:).* exp3 ...
                    -Solution(4+(icr-1)*4,:).*S12(2,:).* exp4;
            else
                S11KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2;
                S22KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2;
                S12KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2;
            end
            S11KW(isnan(S11KW))=0;
            S22KW(isnan(S22KW))=0;
            S12KW(isnan(S12KW))=0;
        end
        
        for ixs=1:nxs
            if nxs~=1
                % integration de la source selon x seulement
                %(gain de temps en calcul)
                % l integration sur la source complete est faite quand
                % nxs==1
                if ~isfield(coordf,'vnx')
                    mxspec=1;
                else
                    if abs(coordf.vnz(ixs))>abs(coordf.vnx(ixs))
                        lsegx=coordf.dr(ixs)*coordf.vnz(ixs);
                        mxspec=sinc(-lsegx*Kpos2/2/pi);
                    else
                        mxspec=1;
                    end
                end
            end
            
            for irx=1:nrecx
                ir = ixr(irx);
                
                if salu(ir)==1
                    %correction spectre
                    U1KW0(1:nk2)            = U1KW(1:nk2).*mxspec;
                    U1KW0(nk2t:-1:(nk2+1))  = U1KW0(2:(nk2t/2)) .*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    U1KW0(1:nk2)            = U1KW0(1:nk2)      .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    UXW1(1,ir,ixs)        	= sum(DK.*U1KW0);
                    
                    %correction spectre
                    U2KW0(1:nk2)            = U2KW(1:nk2).*mxspec;
                    U2KW0(nk2t:-1:(nk2+1))  =-U2KW0(2:(nk2t/2)) .*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    U2KW0(1:nk2)            = U2KW0(1:nk2)      .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    UXW1(2,ir,ixs)          = sum(DK.*U2KW0);
                end
                
                if sals(ir)==1
                    %correction spectre
                    S11KW0(1:nk2)           = S11KW(1:nk2).*mxspec;
                    S11KW0(nk2t:-1:(nk2+1)) = S11KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S11KW0(1:nk2)           = S11KW0(1:nk2)     .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW1(1,1,ir,ixs)        = sum(DK.*S11KW0);
                    
                    %correction spectre
                    S22KW0(1:nk2)           = S22KW(1:nk2).*mxspec;
                    S22KW0(nk2t:-1:(nk2+1)) = S22KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S22KW0(1:nk2)           = S22KW0(1:nk2)     .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW1(2,2,ir,ixs)        = sum(DK.*S22KW0);
                    
                    
                    %correction spectre
                    S12KW0(1:nk2)           = S12KW(1:nk2).*mxspec;
                    S12KW0(nk2t:-1:(nk2+1)) =-S12KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S12KW0(1:nk2)           = S12KW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW1(1,2,ir,ixs)        = sum(DK.*S12KW0);
                end
            end
        end
    end
end


% fint1=0;
% fint2=fij(2);

if fij(2)~=0
    Fsource(is110,:)=Ci(1,1)*G1210 + Ci(1,2)*G2220;                         %sigma11 en 0
    Fsource(is120,:)=Ci(6,6)*(G1220+G2210);                                 %sigma12 en 0
    if ics>1
        Fsource(iu10,:)=G120;                                               %U1 EN 0
        Fsource(iu20,:)=G220;                                               %U2 EN 0
    end
    if (ncapas>1 && ics==1) || (ics<=(ncapas-1))
        Fsource(is11h,:)=Ci(1,1)*G121h + Ci(1,2)*G222h;                     %sigma11 en h
        Fsource(is12h,:)=Ci(6,6)*(G122h+G221h);                             %sigma12 en h
        Fsource(iu1h,:) =G12h;                                              %U1 EN h
        Fsource(iu2h,:) =G22h;                                              %U2 EN h
    end
    Fsource=fij(2)*Fsource*mic;
    
    %Resolution du systeme des conditions aux limites
    for i=1:4*(ncapas-1)+2
        Solution(i,:)   = sum(squeeze(A_DWN(i,:,:)).*Fsource,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des champs au niveau du recepteur %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for irz=1:nrecz
        ixr     = find(izr0==irz);
        nrecx   = length(ixr);
        
        goU=sum(salu(ixr(1:nrecx)));
        goS=sum(sals(ixr(1:nrecx)));
        %Ecriture de la solution en fonction de la position du recepteur
        zricr   = zr0(irz);
        
        icr     = 1;
        while zricr>MAT(icr).h && icr<ncapas
            zricr   = zricr-MAT(icr).h;
            icr     = icr+1;
        end
        
        exp1=exp(-1i.*k1(1,:,icr).* zricr);
        exp2=exp(-1i.*k1(2,:,icr).* zricr);
        exp3=exp( 1i.*k1(1,:,icr).*(zricr-MAT(icr).h));
        exp4=exp( 1i.*k1(2,:,icr).*(zricr-MAT(icr).h));
        
        if goU>0
            %calcul des deplacements diffractes
            if icr<ncapas
                U1KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ...
                    +Solution(3+(icr-1)*4,:).*u1(1,:,icr).* exp3 ...
                    -Solution(4+(icr-1)*4,:).*u1(2,:,icr).* exp4 ;
                U2KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*u2(1,:,icr).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*u2(2,:,icr).* exp4;
            else
                U1KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ;
                U2KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ;
            end
            U1KW(isnan(U1KW))=0;
            U2KW(isnan(U2KW))=0;
        end
        
        if goS>0
            %calcul des contraintes diffractees
            for i=1:2
                S11(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,1).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,2).*u2(i,:,icr));
                S22(i,:)= -1i*(k1(i,:,icr).*MAT(icr).Ci(1,2).*u1(i,:,icr) +k2.*MAT(icr).Ci(1,1).*u2(i,:,icr));
                S12(i,:)= -1i*MAT(icr).Ci(6,6)*( k1(i,:,icr).*u2(i,:,icr) +k2.*u1(i,:,icr));
            end
            if icr<ncapas
                S11KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*S11(1,:).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*S11(2,:).* exp4;
                S22KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2 ...
                    -Solution(3+(icr-1)*4,:).*S22(1,:).* exp3 ...
                    +Solution(4+(icr-1)*4,:).*S22(2,:).* exp4;
                S12KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2 ...
                    +Solution(3+(icr-1)*4,:).*S12(1,:).* exp3 ...
                    -Solution(4+(icr-1)*4,:).*S12(2,:).* exp4;
            else
                S11KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S11(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S11(2,:).* exp2;
                S22KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S22(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S22(2,:).* exp2;
                S12KW(1:nk2) = ...
                    +Solution(1+(icr-1)*4,:).*S12(1,:).* exp1 ...
                    +Solution(2+(icr-1)*4,:).*S12(2,:).* exp2;
            end
            S11KW(isnan(S11KW))=0;
            S22KW(isnan(S22KW))=0;
            S12KW(isnan(S12KW))=0;
        end
        
        
        for ixs=1:nxs
            if nxs~=1
                % integration de la source selon x seulement
                %(gain de temps en calcul)
                % l integration sur la source complete est faite quand
                % nxs==1
                if ~isfield(coordf,'vnx')
                    mxspec=1;
                else
                    if abs(coordf.vnz(ixs))>abs(coordf.vnx(ixs))
                        lsegx=coordf.dr(ixs)*coordf.vnz(ixs);
                        mxspec=sinc(-lsegx*Kpos2/2/pi);
                    else
                        mxspec=1;
                    end
                end
            end
            
            for irx=1:nrecx
                ir = ixr(irx);
                
                if salu(ir)==1
                    %correction spectre
                    U1KW0(1:nk2)            = U1KW(1:nk2).*mxspec;
                    U1KW0(nk2t:-1:(nk2+1))  =-U1KW0(2:(nk2t/2)) .*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    U1KW0(1:nk2)            = U1KW0(1:nk2)      .*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    UXW2(1,ir,ixs)          = sum(DK.*U1KW0);
                    
                    %correction spectre
                    U2KW0(1:nk2)            = U2KW(1:nk2).*mxspec;
                    U2KW0(nk2t:-1:(nk2+1))  = U2KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    U2KW0(1:nk2)            = U2KW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    UXW2(2,ir,ixs)          = sum(DK.*U2KW0);
                end
                
                if sals(ir)==1
                    %correction spectre
                    S11KW0(1:nk2)           = S11KW(1:nk2).*mxspec;
                    S11KW0(nk2t:-1:(nk2+1)) =-S11KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S11KW0(1:nk2)           = S11KW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW2(1,1,ir,ixs)        = sum(DK.*S11KW0);
                    
                    %correction spectre
                    S22KW0(1:nk2)           = S22KW(1:nk2).*mxspec;
                    S22KW0(nk2t:-1:(nk2+1)) =-S22KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S22KW0(1:nk2)           = S22KW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW2(2,2,ir,ixs)        = sum(DK.*S22KW0);
                    
                    %correction spectre
                    S12KW0(1:nk2)           = S12KW(1:nk2).*mxspec;
                    S12KW0(nk2t:-1:(nk2+1)) = S12KW0(2:(nk2t/2)).*exp( 1i*(xr(ir)-xs(ixs))*Kneg2);
                    S12KW0(1:nk2)           = S12KW0(1:nk2).*exp(-1i*(xr(ir)-xs(ixs))*Kpos2);
                    SXW2(1,2,ir,ixs)        = sum(DK.*S12KW0);
                end
            end
        end
    end
end




% ! ojo ici les indices de Gij sont ceux de l'IBEM et
% sont donc inversés par rapport a ceux du DWN
% mais la force finti est dans le repere du DWN
for irz=1:nrecz
    ixr     = find(izr0==irz);
    nrecx   = length(ixr);
    
    %Ecriture de la solution en fonction de la position du recepteur
    zricr   = zr0(irz);
    
    icr     = 1;
    while zricr>MAT(icr).h && icr<ncapas
        zricr   = zricr-MAT(icr).h;
        icr     = icr+1;
    end
    for ixs=1:nxs
        for irx=1:nrecx
            ir = ixr(irx);
            if icr==ics
                xij     = xr(ir)-xs(ixs);
                zij     = zr0(irz)-zs;
                rij     = sqrt(xij.^2+zij.^2);
                ksi 	= para.reg(1).sub(ics).ksi;
                kpi 	= para.reg(1).sub(ics).kpi;
                Ci      = MAT(ics).Ci;
                g(1,1)  = xij./rij;
                g(2,1)  = zij./rij;
                
                if salu(ir)==1
                    if isfield(coordf,'vnx')
                        dr      = coordf.dr(ixs);
                        if rij==0
                            gn(1)=-coordf.vnz(ixs);
                            gn(2)= coordf.vnx(ixs);
                            Gij0 = Greenex_PSV(ksi,kpi,gn,Ci,dr);
                        elseif rij<=1*para.npplo*dr
                            coordf0.z   = zs;
                            coordf0.x   = xs(ixs);
                            coordf0.vnx = coordf.vnx(ixs);
                            coordf0.vnz = coordf.vnz(ixs);
                            coordf0.dr  = dr;
                            Gij0 = Gij_PSV_r_small(coordf0,xr(ir),zr0(irz),1,ksi,kpi,para.gaussian,Ci);
                        else
                            Gij0 = Gij_PSV(ksi,kpi,rij,g,Ci,1);
                        end
                    else
                        Gij0 = Gij_PSV(ksi,kpi,rij,g,Ci,1);
                    end
                    UXW1(1,ir,ixs) = UXW1(1,ir,ixs) + fij(1)*Gij0(2,2);
                    UXW1(2,ir,ixs) = UXW1(2,ir,ixs) + fij(1)*Gij0(1,2);
                    UXW2(1,ir,ixs) = UXW2(1,ir,ixs) + fij(2)*Gij0(1,2);
                    UXW2(2,ir,ixs) = UXW2(2,ir,ixs) + fij(2)*Gij0(1,1);
                    
                end
                
                if sals(ir)==1
                    if isfield(coordf,'vnx')
                        dr      = coordf.dr(ixs);
                        if rij<=1*para.npplo*dr && rij~=0
                            coordf0.z   = zs;
                            coordf0.x   = xs(ixs);
                            coordf0.vnx = coordf.vnx(ixs);
                            coordf0.vnz = coordf.vnz(ixs);
                            coordf0.dr  = dr;
                            
                            [S2,S1]=S_PSV_r_small_2(coordf0,xr(ir),zr0(irz),ksi,kpi,Ci,fijIBEM,para.gaussian);
                            
                            SXW1(1,1,ir,ixs) = SXW1(1,1,ir,ixs) + fijIBEM(2)*S1(2,2);%ti=sij.nj
                            SXW1(2,2,ir,ixs) = SXW1(2,2,ir,ixs) + fijIBEM(2)*S1(1,1);
                            SXW1(1,2,ir,ixs) = SXW1(1,2,ir,ixs) + fijIBEM(2)*S1(1,2);
                            
                            SXW2(1,1,ir,ixs) = SXW2(1,1,ir,ixs) + fijIBEM(1)*S2(2,2);%ti=sij.nj
                            SXW2(2,2,ir,ixs) = SXW2(2,2,ir,ixs) + fijIBEM(1)*S2(1,1);
                            SXW2(1,2,ir,ixs) = SXW2(1,2,ir,ixs) + fijIBEM(1)*S2(1,2);
                        elseif rij~=0
                            %dans le repere IBEM
                            [S2,S1]=S_PSV(rij,g,ksi,kpi,Ci,fijIBEM); %a ameliorer, on passe 2 fois pour chacun des fij et on peut le faire dehors
                            SXW1(1,1,ir,ixs) = SXW1(1,1,ir,ixs) + fijIBEM(2)*S1(2,2);%ti=sij.nj
                            SXW1(2,2,ir,ixs) = SXW1(2,2,ir,ixs) + fijIBEM(2)*S1(1,1);
                            SXW1(1,2,ir,ixs) = SXW1(1,2,ir,ixs) + fijIBEM(2)*S1(1,2);
                            
                            SXW2(1,1,ir,ixs) = SXW2(1,1,ir,ixs) + fijIBEM(1)*S2(2,2);%ti=sij.nj
                            SXW2(2,2,ir,ixs) = SXW2(2,2,ir,ixs) + fijIBEM(1)*S2(1,1);
                            SXW2(1,2,ir,ixs) = SXW2(1,2,ir,ixs) + fijIBEM(1)*S2(1,2);
                        end
                    else
                        %dans le repere IBEM
                        [S2,S1]=S_PSV(rij,g,ksi,kpi,Ci,fijIBEM); %a ameliorer, on passe 2 fois pour chacun des fij et on peut le faire dehors
                        % avec chgmt d indice de Gij, fint ok
                        
                        SXW1(1,1,ir,ixs) = SXW1(1,1,ir,ixs) + fijIBEM(2)*S1(2,2);%ti=sij.nj
                        SXW1(2,2,ir,ixs) = SXW1(2,2,ir,ixs) + fijIBEM(2)*S1(1,1);
                        SXW1(1,2,ir,ixs) = SXW1(1,2,ir,ixs) + fijIBEM(2)*S1(1,2);
                        
                        SXW2(1,1,ir,ixs) = SXW2(1,1,ir,ixs) + fijIBEM(1)*S2(2,2);%ti=sij.nj
                        SXW2(2,2,ir,ixs) = SXW2(2,2,ir,ixs) + fijIBEM(1)*S2(1,1);
                        SXW2(1,2,ir,ixs) = SXW2(1,2,ir,ixs) + fijIBEM(1)*S2(1,2);
                    end
                end
            end
        end
    end
end
%cuidado 1 y 2 son intercambiados entre IBEM y DWN
UXW1(1:2,:,:)= UXW1(2:-1:1,:,:);
SXW0         = SXW1;
SXW1(1,1,:,:)= SXW0(2,2,:,:);
SXW1(2,2,:,:)= SXW0(1,1,:,:);

UXW2(1:2,:,:)= UXW2(2:-1:1,:,:);
SXW0         = SXW2;
SXW2(1,1,:,:)= SXW0(2,2,:,:);
SXW2(2,2,:,:)= SXW0(1,1,:,:);