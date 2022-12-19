function [UXW,SXW]=calcul_US_DWN_SH_Ncapas_HS_incOP(para,xr,zr0,izr0,salu,sals,coordf,kxk,kzsigno)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)
%           SXW = [SxyKW,SzyKW](xr,zr,w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% attention, pour calculer les differentes composantes
% il est nécessaire de subdiviser le calcul en fonction
% des parties paire et impaire (k2)
% cf G11, G22 paires et G12 impaire
% ce qui peut être fait en traitant séparement fint1 et fint2
nk2         = length(kxk);

ncapas      = para.nsubmed;
A_DWN       = zeros(2*(ncapas-1)+1,2*(ncapas-1)+1,nk2);
nrec        = length(xr);
nrecz       = length(zr0);

UyKW        = 0;
SzyKW       = 0;
SxyKW       = 0;

xs          = coordf.xs;
nxs       	= length(xs);
zs          = coordf.zs;

UXW         = zeros(nrec,nxs);
SXW         = zeros(2,nrec,nxs);

Fsource     = zeros(2*(ncapas-1)+1,nk2);

% identification des composantes du nombre d onde incident
% l incidence provient du semi-espace
ics = ncapas;
ksi	= para.reg(1).sub(ncapas).ksi;
% ksi	= para.reg(1).sub(1).ksi;
k2  = kxk*ksi;
k1f	= kzsigno.*sqrt(ksi^2-k2.^2); %kz = + o -, la 2 soluciones son correctas


%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
k1 = zeros(ncapas,nk2);
for ic=1:ncapas
    ksi     = para.reg(1).sub(ic).ksi;
    k1(ic,:)  = sqrt(ksi^2-k2.^2); %k1S
end
k1 =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));
% DWN.k10=k1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%les ondes en exp(-1i.*k1(ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*k1(ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %*********************************************************%
    %   calcul de la contrainte issue de la solution homogene %
    %*********************************************************%
    MAT = para.reg(1).sub(ic);
    S12 =-1i*k1(ic,:)*MAT.Ci(6,6); %sigma_zy sigma12
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    %SIGMA12 EN 0
    A_DWN(1,1,:)= mic*S12; %.*exp(-1i.*k1(ic).*0);                        %Amplitude trans sens +
    if ncapas>1
        exp2=exp(-1i.*k1(ic,:).*MAT.h);
        exp4=exp( 1i.*k1(ic,:).*(0-MAT.h));
        
        A_DWN(1,2,:)=-mic*S12.*exp4; %.*exp( 1i.*k1(ic).*0);                    %Amplitude trans sens -
        
        %SIGMA12 EN h(ic)
        A_DWN(2,1,:)= mic* S12.*exp2;                                        %Amplitude trans sens +
        A_DWN(2,2,:)=-mic* S12;                                        %Amplitude trans sens -
        
        %U2 EN h(ic)
        A_DWN(3,1,:)= mic.*exp2;                                          %Amplitude trans sens +
        A_DWN(3,2,:)= mic;                                          %Amplitude trans sens -
    end
end

for ic=2:(ncapas-1)
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    MAT = para.reg(1).sub(ic);
    S12 =-1i*k1(ic,:)*MAT.Ci(6,6); %sigma_12
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    exp2= exp(-1i.*k1(ic,:).*MAT.h);
    exp4= exp( 1i.*k1(ic,:).*(0-MAT.h));
    mic	= (-1)^ic;
    
    %SIGMA12 EN 0
    A_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*S12; %.*exp(-1i.*k1(ic).*0);	%Amplitude trans sens +
    A_DWN((ic-2)*2+1+1,(ic-1)*2+2,:)=-mic*S12.*exp4; %.*exp( 1i.*k1(ic).*0);%Amplitude trans sens -
    
    %U2 EN 0
    A_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic; %.*exp(-1i.*k1(ic).*0);        %Amplitude trans sens +
    A_DWN((ic-2)*2+1+2,(ic-1)*2+2,:)= mic.*exp4; %.*exp( 1i.*k1(ic).*0);	%Amplitude trans sens -
    
    %SIGMA12 EN h(ic)
    A_DWN((ic-2)*2+1+3,(ic-1)*2+1,:)= mic*S12.*exp2;                        %Amplitude trans sens +
    A_DWN((ic-2)*2+1+3,(ic-1)*2+2,:)=-mic*S12;                              %Amplitude trans sens -
    
    %U2 EN h(ic)
    A_DWN((ic-2)*2+1+4,(ic-1)*2+1,:)= mic.*exp2;                            %Amplitude trans sens +
    A_DWN((ic-2)*2+1+4,(ic-1)*2+2,:)= mic;                                  %Amplitude trans sens -
end

if ncapas>1
    for ic=ncapas
        %***********************************************************%
        %   calcul de la contrainte issue de la solution homogene   %
        %***********************************************************%
        MAT = para.reg(1).sub(ic);
        S12 =-1i*k1(ic,:)*MAT.Ci(6,6); %sigma_12
        
        mic	= (-1)^ic;
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        %SIGMA12 EN 0
        A_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*S12; %.*exp(-1i.*k1(ic).*0);		%Amplitude trans sens +
        
        %U2 EN 0
        A_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic; %.*exp(-1i.*k1(ic).*0);	%Amplitude trans sens +
    end
end




MAT = para.reg(1).sub;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul du vecteur source et resolution des coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% le calcul du vecteur source se fait en utilisant
% la decomposition en ondes planes pour l element j
if ics>1
    htop=0;
    for inch=1:ics-1
        htop=htop+MAT(inch).h;
    end
else
    htop=0;
end

%%%%%%%%%%%%%%%%%
% interface sup %
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des deplacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
zt      = htop-zs;%x=z dans l epaisseur
exp2t   = exp(1i.*k1f.*zt);
G220    = exp2t;%pair en k2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des derivees utiles au calcul des contraintes normales %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivee p/x1 en X (x1=z)
%pair
G22z0=1i*k1f.*exp2t;

%indice de los terminos en el vector fuente
if ics==1
    is120=1;%SIGMA12 EN 0
elseif ics>1
    is120=(ics-2)*2+1+1;%SIGMA12 EN 0
    iu20=(ics-2)*2+1+2;%U2 EN 0
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci  = MAT(ics).Ci;
Fsource(is120,:)=Ci(6,6)*G22z0;                           %sigma12 en 0
if ics>1
    Fsource(iu20,:)=G220;%U2 EN 0
end
mic     = -(-1)^ics;
Fsource = Fsource*mic;

%Resolution du systeme des conditions aux limites
Solution=zeros(2*(ncapas-1)+1,nk2);
for i=1:nk2
    Solution(:,i)   = A_DWN(:,:,i)\Fsource(:,i);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % normalisation en energie %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E = 0;
% % for icr=1:ncapas-1
% %     %negligeable
% % end
% icr = ncapas;
% E   = E + ( (1+ abs(Solution(1+(icr-1)*2))^2)/2 -0*imag(Solution(1+(icr-1)*2)) )./k1f;
% E = sqrt(real(E));
% 
% for icr=1:ncapas-1
%     Solution(1+(icr-1)*2)=Solution(1+(icr-1)*2)/E;
%     Solution(2+(icr-1)*2)=Solution(2+(icr-1)*2)/E;
% end
% icr = ncapas;
% Solution(1+(icr-1)*2)=Solution(1+(icr-1)*2)/E;
E=1;

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
    
    exp1=exp(-1i.*k1(icr,:).*zricr).';
    exp3=exp( 1i.*k1(icr,:).*(zricr-MAT(icr).h)).';
    
    if goU>0
        %calcul des deplacements diffractes
        if icr<ncapas
            UyKW = ...
                +Solution(1+(icr-1)*2,:).* exp1.' ...
                +Solution(2+(icr-1)*2,:).* exp3.';
        else
            UyKW = Solution(1+(icr-1)*2,:).* exp1.';
        end
        UyKW(isnan(UyKW))=0;
    end
    
    if goS>0
        %calcul des contraintes diffractees
        %pair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k1(icr,:).';
        if icr<ncapas
            SzyKW = ...
                +Solution(1+(icr-1)*2,:).*C0.* exp1.' ...
                -Solution(2+(icr-1)*2,:).*C0.* exp3.';
        else
            SzyKW =Solution(1+(icr-1)*2).*C0.* exp1.';
        end
        SzyKW(isnan(SzyKW))=0;
        
        %impair en k2
        C0= -1i*MAT(icr).Ci(6,6)*k2;
        if icr<ncapas
            SxyKW = ...
                +Solution(1+(icr-1)*2,:).*C0.* exp1.' ...
                +Solution(2+(icr-1)*2,:).*C0.* exp3.';
        else
            SxyKW =Solution(1+(icr-1)*2).*C0.* exp1.';
        end
        SxyKW(isnan(SxyKW))=0;
    end
    
    %rajouts des champs incidents
    if icr==ics
        zrs=zr0(irz)-zs;
        exp2t= exp(1i.*k1f.*zrs)/E;
        if goU>0
            G22r = exp2t;
            
            %calcul des deplacements incidents
            UyKW = UyKW + G22r;
        end
        
        if goS>0
            %%%derivee p/x1 en X (x1=z)
            %pair
            G22z=1i*k1f.*exp2t.';
            
            % derivee p/x2 en X (x2=x)
            %impair
            G22x=-1i*k2.*exp2t.';
            
            %calcul des tractions t1j
            SzyKW = SzyKW + MAT(icr).Ci(6,6)*G22z;
            SxyKW = SxyKW + MAT(icr).Ci(6,6)*G22x;
        end
    end
    
    if min(xs==xs(1))==1 %todos son iguales
        for irx=1:nrecx
            ir  = ixr(irx);
            expx= exp(-1i*(xr(ir)-xs).*k2);
            if salu(ir)==1
                UXW(ir,:) = UyKW.*expx;
            end
            
            if sals(ir)==1
                SXW(1,ir,:)=SxyKW.*expx;
                SXW(2,ir,:)=SzyKW.*expx;
            end
        end
    else
        for ixs=1:nxs
            for irx=1:nrecx
                ir  = ixr(irx);
                expx= exp(-1i*(xr(ir)-xs(ixs))*k2(ixs));
                if salu(ir)==1
                    UXW(ir,ixs) = UyKW.*expx;
                end
                
                if sals(ir)==1
                    SXW(1,ir,ixs)=SxyKW.*expx;
                    SXW(2,ir,ixs)=SzyKW.*expx;
                end
            end
        end
    end
end