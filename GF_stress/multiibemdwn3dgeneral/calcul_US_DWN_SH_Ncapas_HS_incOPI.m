function [UXW,SXW]=calcul_US_DWN_SH_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno)
% caso de incidencia de ondas inhomogenea en multi estratos
% en este caso se supone por lo menos 1 estrato sobre un semi espacio, el
% caso del semi-espacio no tiene ondas inhomogenea
% Se impone una incidencia de amplitud unitaria
% la onda incidente se impone en el estrato de uno de los receptores para
% evitar problemas numericos
% se busca las amplitudes de las ondas en los estratos tomando en cuenta la
% incidencia de la onda que se impuso
% asi el sistema se reduce de una dimension con respeto al caso normal
% se quita la linea s12 donde esta la fuente, y la amplitud de la onda inc
% se normaliza la energia de cada modo de manera que la integral con
% respecto a la profundidad es unitaria


ncapas      = para.nsubmed;
nrec        = length(xr);
nrecz       = length(zr0);


% zs          = coordf.zs;

%identificacion del estrato del receptor de minima profundidad
zrm         = min(zr0);
ics         = 1;
MAT         = para.reg(1).sub;
hc          = MAT(1).h;
while zrm>hc && ics<ncapas
    ics     = ics+1;
    hc      = hc+MAT(ics).h;
end
% se considera que la onda incidente proviene del estrato del receptor de
% minima profundidad
xs          = mean(xr);
nxs       	= length(xs); %### a desarollar para varias incidencias

%inicialisaciones
UyKW        = 0;
SzyKW       = 0;
SxyKW       = 0;
UXW         = zeros(nrec,nxs);
SXW         = zeros(2,nrec,nxs);
Fsource     = zeros(2*(ncapas-1),1);
A_DWN       = zeros(2*(ncapas-1),2*(ncapas-1));

% identification des composantes du nombre d onde incident
% l incidence provient du semi-espace (par convention)
ksi         = para.reg(1).sub(ncapas).ksi;
k2          = kxk*ksi;

%**********************************************************************%
%  calcul des nombres d'onde et polarisations de la solution homogene  %
%**********************************************************************%
k1 = zeros(ncapas,1);
for ic=1:ncapas
    ksi     = para.reg(1).sub(ic).ksi;
    k1(ic)  = sqrt(ksi^2-k2.^2); %k1S
end
k1 =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));

%************************************************%
% calcul de la matriz de condicion de interfaces %
%************************************************%
%les ondes en exp(-1i.*k1(ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*k1(ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %*********************************************************%
    %   calcul de la contrainte issue de la solution homogene %
    %*********************************************************%
    mic                 = (-1)^ic;%convention de signe entre les equations AX=B
    MAT                 = para.reg(1).sub(ic);
    S12                 =-1i*k1(ic)*MAT.Ci(6,6); %sigma_zy == sigma12
    %translation d indice pour assurer la construction de la matrice sans
    %prendre en compte l incidence imposee
    itr                 = real(ics<=ic);

    %SIGMA12 EN 0
    if ics~=ic
        exp2            = exp(-1i.*k1(ic).*MAT.h);
        exp4          	= exp( 1i.*k1(ic).*(0-MAT.h));
        A_DWN(1,1)      = mic*S12; %.*exp(-1i.*k1(ic).*0);                  %Amplitude trans sens +
        A_DWN(1,2-itr)	=-mic*S12.*exp4; %.*exp( 1i.*k1(ic).*0);            %Amplitude trans sens -
    end
    
    %SIGMA12 EN h(ic)
    if ics~=ic+1
        if ics~=ic
            A_DWN(2,1)      = mic* S12.*exp2;                                   %Amplitude trans sens +
        end
        A_DWN(2-itr,2-itr)  =-mic* S12;                                         %Amplitude trans sens -
        itr1=itr;
    else
        itr1=itr+1;
    end
    
    %U2 EN h(ic)
    if ics~=ic
        A_DWN(3-itr1,1)  = mic.*exp2;                                        %Amplitude trans sens +
    end
    A_DWN(3-itr1,2-itr)  = mic;                                          %Amplitude trans sens -
end

for ic=2:(ncapas-1)
    mic	= (-1)^ic;
    itr = real(ics<=ic);
    MAT = para.reg(1).sub(ic);
    S12 =-1i*k1(ic)*MAT.Ci(6,6); %sigma_12
    exp2= exp(-1i.*k1(ic).*MAT.h);
    exp4= exp( 1i.*k1(ic).*(0-MAT.h));
    
    %SIGMA12 EN 0
    if ics~=ic
        A_DWN((ic-2)*2+1+1-itr,(ic-1)*2+1-itr)= mic*S12; %.*exp(-1i.*k1(ic).*0);	%Amplitude trans sens +
        A_DWN((ic-2)*2+1+1-itr,(ic-1)*2+2-itr)=-mic*S12.*exp4; %.*exp( 1i.*k1(ic).*0);%Amplitude trans sens -
    end
    
    %U2 EN 0
    if ics~=ic
        A_DWN((ic-2)*2+1+2-itr,(ic-1)*2+1-itr)= mic; %.*exp(-1i.*k1(ic).*0);        %Amplitude trans sens +
    end
    A_DWN((ic-2)*2+1+2-itr,(ic-1)*2+2-itr)= mic.*exp4; %.*exp( 1i.*k1(ic).*0);	%Amplitude trans sens -
    
    %SIGMA12 EN h(ic)
    if ics~=ic+1
        if ics~=ic
            A_DWN((ic-2)*2+1+3-itr,(ic-1)*2+1-itr)= mic*S12.*exp2;                        %Amplitude trans sens +
        end
        A_DWN((ic-2)*2+1+3-itr,(ic-1)*2+2-itr)=-mic*S12;                              %Amplitude trans sens -
         itr1=itr;
    else
        itr1=itr+1;
    end
    
    %U2 EN h(ic)
    if ics~=ic
        A_DWN((ic-2)*2+1+4-itr1,(ic-1)*2+1-itr)= mic.*exp2;                            %Amplitude trans sens +
    end
    A_DWN((ic-2)*2+1+4-itr1,(ic-1)*2+2-itr)= mic;                                  %Amplitude trans sens -
end

if ics~=ncapas
    ic  = ncapas;
    itr = real(ics<=ic);

    mic	= (-1)^ic;
    itr = real(ics<=ic);
    S12 =-1i*k1(ic)*para.reg(1).sub(ic).Ci(6,6); %sigma_12
    
    %SIGMA12 EN 0
    A_DWN((ic-2)*2+1+1-itr,(ic-1)*2+1-itr)= mic*S12; %.*exp(-1i.*k1(ic).*0);		%Amplitude trans sens +
    %U2 EN 0
    A_DWN((ic-2)*2+1+2-itr,(ic-1)*2+1-itr)= mic; %.*exp(-1i.*k1(ic).*0);	%Amplitude trans sens +
end

MAT = para.reg(1).sub;

%**************************%
% calcul du vecteur source %
%**************************%

%indice de los terminos en el vector fuente
if ics==1
    is12h=1;    %SIGMA12 EN h
    iu2h =2;    %U2 EN h
elseif ics==ncapas
    iu20 =(ics-2)*2+1+1;%U2 EN 0
else %1<ics<ncapas
    iu20 =(ics-2)*2+1+1;%U2 EN 0
    is12h=(ics-2)*2+1+2;%SIGMA12 EN h
    iu2h =(ics-2)*2+1+3;%U2 EN h
end

% calcul du terme source %
k1f    	= k1(ics);%k1S
Ci      = MAT(ics).Ci;
u2inc   = 1;
S12     =-1i*k1f*Ci(6,6);
exp2t   = exp(-1i.*k1f.*MAT(ics).h);

%Fsource(is120,:)=S12;                    	%sigma12 en 0
if ics>1
    Fsource(iu20 ,:)= u2inc;                %u2 EN 0
end                 
if ics<ncapas
    Fsource(is12h,:)= S12  .*exp2t;      	%sigma12 en h
    Fsource(iu2h ,:)= u2inc.*exp2t;      	%U2 EN h
end

mic     = -(-1)^ics;
Fsource = Fsource*mic;

%*****************************%
% resolution des coefficients %
%*****************************%

%Resolution du systeme des conditions aux limites
% Solution   = A_DWN\Fsource;
Solution    = A_DWN\Fsource;
% Solution    = pinv(A_DWN)*Fsource;

%on reinsere la composante incidente ds la solution
%et on supprime plus tard la contribution du champ incident
n0 = 2*(para.nsubmed-1);
nk2=1;
for i=1:nk2
    Solution((((ics-1)*2+1):n0)+1,i) = Solution(((ics-1)*2+1):n0,i);
end
Solution((ics-1)*2+1,i)   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalisation en energie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%l origine de l'onde ds le SE doit etre au niveau de l interface
%SE/derniere couche pour que la normalisation soit correcte
E = 0;
for icr=1:ncapas-1
    %     k1i = imag(k1);
    %     exp1= exp(-1i*2*k1(icr).*MAT(icr).h).';
    %     exp3= exp( 1i*2*k1(icr).*(0-MAT(icr).h)).';
    S1  = Solution(1+(icr-1)*2);
    S2  = Solution(2+(icr-1)*2);
    h   = MAT(icr).h;
    kz  = k1(icr);
    
    %k1 real
    %     E   = E + (abs(S1)^2 + abs(S2)^2 )*h + (S1*conj(S2) + S2*conj(S1))*sin(kz*h)/kz;
    
    %k1 complex
    E   = E + kz*1i*(...
        (S1*conj(S2)+S2*conj(S1))*(exp(-h*kz*1i) - exp(h*abs(kz)^2*1i/kz))/(kz^2 + abs(kz)^2) + ...
        (S1*conj(S1)+S2*conj(S2))*(exp(h*1i*((abs(kz)^2)/kz-kz)) - 1)     /(kz^2 - abs(kz)^2));
end
% icr = ncapas;
k1f = imag(k1(ncapas));
E   = E - 1/(2*k1f);
E   = real(sqrt(E));

%normalisation des amplitudes
Solution= Solution/E;


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
        zricr  = zricr-MAT(icr).h;
        icr = icr+1;
    end
    
    exp1=exp(-1i.*k1(icr).*zricr).';
    exp3=exp( 1i.*k1(icr).*(zricr-MAT(icr).h)).';
    
    if goU>0
        %calcul des deplacements diffractes
        if icr<ncapas
            UyKW = ...
                +Solution(1+(icr-1)*2).* exp1 ...
                +Solution(2+(icr-1)*2).* exp3;
        else
            UyKW = Solution(1+(icr-1)*2).* exp1 ;
        end
        UyKW(isnan(UyKW))=0;
    end
    
    if goS>0
        %calcul des contraintes diffractees
        %pair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k1(icr).';
        if icr<ncapas
            SzyKW = ...
                +Solution(1+(icr-1)*2).*C0.* exp1 ...
                -Solution(2+(icr-1)*2).*C0.* exp3;
        else
            SzyKW = Solution(1+(icr-1)*2).*C0.* exp1;
        end
        SzyKW(isnan(SzyKW))=0;
        
        %impair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k2;
        if icr<ncapas
            SxyKW = ...
                +Solution(1+(icr-1)*2).*C0.* exp1 ...
                +Solution(2+(icr-1)*2).*C0.* exp3;
        else
            SxyKW = Solution(1+(icr-1)*2).*C0.* exp1;
        end
        SxyKW(isnan(SxyKW))=0;
    end
    
    for ixs=1:nxs
        for irx=1:nrecx
            ir  = ixr(irx);
            expx= exp(-1i*(xr(ir)-xs(ixs))*k2);
            if salu(ir)==1
                UXW(ir,ixs) = UyKW*expx;
            end
            
            if sals(ir)==1
                SXW(1,ir,ixs)=SxyKW*expx;
                SXW(2,ir,ixs)=SzyKW*expx;
            end
        end
    end
end