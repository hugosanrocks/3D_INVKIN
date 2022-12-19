function [UXW,SXW]=calcul_US_DWN_PSV_Ncapas_HS_incOPIP(para,xr,zr0,izr0,salu,sals,coordf,kxk)
% METHODE : @ recuperacion de k1, u1, u2 (sol homogenea)
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)
% se quita la onda reflejada en el semi-espacio
% se busca las amplitudes de las ondas en los estratos con una mimimizacion
% de las ecuaciones en las interfases (ya que hay una ecuacion mas que de
% incognitas)
% se normaliza la energia de cada modo de manera que la integral con
% respecto a la profundidad es unitaria
% se multiplica por su factor de particion


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% attention, pour calculer les differentes composantes
% il est nécessaire de subdiviser le calcul en fonction
% des parties paire et impaire (k2)
% cf G11, G22 paires et G12 impaire
% ce qui peut être fait en traitant séparement fint1 et fint2
ncapas      = para.nsubmed;
nk2         = length(kxk);
nrec        = length(xr);
nrecz       = length(zr0);

xs          = coordf.xs;
nxs       	= length(xs);
zs          = coordf.zs;

UXW         = zeros(2,nrec,nxs);
SXW         = zeros(2,2,nrec,nxs);


S22         = zeros(2,1);
Fsource     = zeros(4*(ncapas-1)+1,nk2);
% Solution    = zeros(4*(ncapas-1)+2,1);

% identification des composantes du nombre d onde incident
% l incidence provient du semi-espace
ics         = ncapas;
ksi         = para.reg(1).sub(ncapas).ksi;
k2          = kxk*ksi;

%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
k1              = zeros(2,nk2,ncapas);
for ic=1:ncapas
    ksi         = para.reg(1).sub(ic).ksi;
    kpi         = para.reg(1).sub(ic).kpi;
    k1(1,:,ic)  = sqrt(ksi^2-k2.^2); %k1S
    k1(2,:,ic)  = sqrt(kpi^2-k2.^2); %k1P
end
k1              =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));

u1              = zeros(2,nk2,ncapas);
u2              = zeros(2,nk2,ncapas);

S11             = zeros(2,nk2);
S12             = zeros(2,nk2);

A_DWN           = zeros(4*(para.nsubmed-1)+1,4*(para.nsubmed-1)+1,nk2);

for ic=1:ncapas
    u1(1,:,ic)  =  k2;          % S
    u2(1,:,ic)  = -k1(1,:,ic);  % S
    u1(2,:,ic)  =  k1(2,:,ic);  % P
    u2(2,:,ic)  =  k2;          % P
end

un      = sqrt(abs(u1).^2+abs(u2).^2);
u1      = u1./un;
u2      = u2./un;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%les ondes en exp(-1i.*k1(:,ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*k1(:,ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    MAT=para.reg(1).sub(ic);
    
    for i=1:2
        S11(i,:)= -1i*(k1(i,:,ic).*MAT.Ci(1,1).*u1(i,:,ic) +k2.*MAT.Ci(1,2).*u2(i,:,ic));   %sigma_11 %Amplitude longi & transverse sens z + (profondeur)
        S12(i,:)= -1i*MAT.Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));                %sigma_12 %Amplitude longi & transverse sens z +
    end
    % S11(1,-)=-S11(1,+) pour k1-=-k1+ & S
    % S12(1,-)= S12(1,+)
    % S11(2,-)= S11(2,+) pour k1-=-k1+ & P
    % S12(2,-)=-S12(2,+)
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    A_DWN(1,1,:)=mic*S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);                 %Amplitude trans sens +
    A_DWN(1,2,:)=mic*S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);                 %Amplitude longi sens +
    %SIGMA12 EN 0
    A_DWN(2,1,:)=mic*S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);                 %Amplitude trans sens +
    A_DWN(2,2,:)=mic*S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);                 %Amplitude longi sens +
    
    exp1=exp(-1i.*k1(1,:,ic).*MAT.h);
    exp2=exp(-1i.*k1(2,:,ic).*MAT.h);
    exp3=exp( 1i.*k1(1,:,ic).*(0-MAT.h));
    exp4=exp( 1i.*k1(2,:,ic).*(0-MAT.h));
    
    %SIGMA11 EN 0
    A_DWN(1,3,:)=-mic*S11(1,:).*exp3;                                   %Amplitude trans sens -
    A_DWN(1,4,:)= mic*S11(2,:).*exp4;                                   %Amplitude longi sens -
    %SIGMA12 EN 0
    A_DWN(2,3,:)= mic*S12(1,:).*exp3;                                   %Amplitude trans sens -
    A_DWN(2,4,:)=-mic*S12(2,:).*exp4;                                   %Amplitude longi sens -
    
    %SIGMA11 EN h(ic)
    A_DWN(3,1,:)=mic* S11(1,:).*exp1;                                   %Amplitude trans sens +
    A_DWN(3,2,:)=mic* S11(2,:).*exp2;                                   %Amplitude longi sens +
    A_DWN(3,3,:)=mic*-S11(1,:);                                         %Amplitude trans sens -
    A_DWN(3,4,:)=mic* S11(2,:);                                         %Amplitude longi sens -
    
    %SIGMA12 EN h(ic)
    A_DWN(4,1,:)=mic* S12(1,:).*exp1;                                   %Amplitude trans sens +
    A_DWN(4,2,:)=mic* S12(2,:).*exp2;                                   %Amplitude longi sens +
    A_DWN(4,3,:)=mic* S12(1,:);                                         %Amplitude trans sens -
    A_DWN(4,4,:)=mic*-S12(2,:);                                         %Amplitude longi sens -
    
    %U1 EN h(ic)
    A_DWN(5,1,:)=mic* u1(1,:,ic).*exp1;                                 %Amplitude trans sens +
    A_DWN(5,2,:)=mic* u1(2,:,ic).*exp2;                                 %Amplitude longi sens +
    A_DWN(5,3,:)=mic* u1(1,:,ic);                                       %Amplitude trans sens -
    A_DWN(5,4,:)=mic*-u1(2,:,ic);                                       %Amplitude longi sens -
    
    if ncapas>2
        %U2 EN h(ic)
        A_DWN(6,1,:)=mic* u2(1,:,ic).*exp1;                                 %Amplitude trans sens +
        A_DWN(6,2,:)=mic* u2(2,:,ic).*exp2;                                 %Amplitude longi sens +
        A_DWN(6,3,:)=mic*-u2(1,:,ic);                                       %Amplitude trans sens -
        A_DWN(6,4,:)=mic* u2(2,:,ic);                                       %Amplitude longi sens -
    end
end

for ic=2:(ncapas-1)
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    MAT=para.reg(1).sub(ic);
    
    for i=1:2
        S11(i,:)= -1i*(k1(i,:,ic).*MAT.Ci(1,1).*u1(i,:,ic) +k2.*MAT.Ci(1,2).*u2(i,:,ic));
        S12(i,:)= -1i*MAT.Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
    end
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    exp1= exp(-1i.*k1(1,:,ic).*MAT.h);
    exp2= exp(-1i.*k1(2,:,ic).*MAT.h);
    exp3= exp( 1i.*k1(1,:,ic).*(0-MAT.h));
    exp4= exp( 1i.*k1(2,:,ic).*(0-MAT.h));
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)=mic* S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    A_DWN((ic-2)*4+2+1,(ic-1)*4+2,:)=mic* S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    A_DWN((ic-2)*4+2+1,(ic-1)*4+3,:)=mic*-S11(1,:).*exp3;                           %Amplitude trans sens -
    A_DWN((ic-2)*4+2+1,(ic-1)*4+4,:)=mic* S11(2,:).*exp4;                           %Amplitude longi sens -
    
    %SIGMA12 EN 0
    A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    A_DWN((ic-2)*4+2+2,(ic-1)*4+2,:)=mic* S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    A_DWN((ic-2)*4+2+2,(ic-1)*4+3,:)=mic* S12(1,:).*exp3;                           %Amplitude trans sens -
    A_DWN((ic-2)*4+2+2,(ic-1)*4+4,:)=mic*-S12(2,:).*exp4;                           %Amplitude longi sens -
    
    %U1 EN 0
    A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    A_DWN((ic-2)*4+2+3,(ic-1)*4+2,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    A_DWN((ic-2)*4+2+3,(ic-1)*4+3,:)=mic* u1(1,:,ic).*exp3;                       	%Amplitude trans sens -
    A_DWN((ic-2)*4+2+3,(ic-1)*4+4,:)=mic*-u1(2,:,ic).*exp4;                         %Amplitude longi sens -
    
    %U2 EN 0
    A_DWN((ic-2)*4+2+4,(ic-1)*4+1,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    A_DWN((ic-2)*4+2+4,(ic-1)*4+2,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    A_DWN((ic-2)*4+2+4,(ic-1)*4+3,:)=mic*-u2(1,:,ic).*exp3;                         %Amplitude trans sens -
    A_DWN((ic-2)*4+2+4,(ic-1)*4+4,:)=mic* u2(2,:,ic).*exp4;                         %Amplitude longi sens -
    
    %SIGMA11 EN h(ic)
    A_DWN((ic-2)*4+2+5,(ic-1)*4+1,:)=mic* S11(1,:).*exp1;                   %Amplitude trans sens +
    A_DWN((ic-2)*4+2+5,(ic-1)*4+2,:)=mic* S11(2,:).*exp2;                   %Amplitude longi sens +
    A_DWN((ic-2)*4+2+5,(ic-1)*4+3,:)=mic*-S11(1,:);                         %Amplitude trans sens -
    A_DWN((ic-2)*4+2+5,(ic-1)*4+4,:)=mic* S11(2,:);                         %Amplitude longi sens -
    
    %SIGMA12 EN h(ic)
    A_DWN((ic-2)*4+2+6,(ic-1)*4+1,:)=mic* S12(1,:).*exp1;                   %Amplitude trans sens +
    A_DWN((ic-2)*4+2+6,(ic-1)*4+2,:)=mic* S12(2,:).*exp2;                   %Amplitude longi sens +
    A_DWN((ic-2)*4+2+6,(ic-1)*4+3,:)=mic* S12(1,:);                         %Amplitude trans sens -
    A_DWN((ic-2)*4+2+6,(ic-1)*4+4,:)=mic*-S12(2,:);                         %Amplitude longi sens -
    
    %U1 EN h(ic)
    A_DWN((ic-2)*4+2+7,(ic-1)*4+1,:)=mic* u1(1,:,ic).*exp1;                 %Amplitude trans sens +
    A_DWN((ic-2)*4+2+7,(ic-1)*4+2,:)=mic* u1(2,:,ic).*exp2;                 %Amplitude longi sens +
    A_DWN((ic-2)*4+2+7,(ic-1)*4+3,:)=mic* u1(1,:,ic);                       %Amplitude trans sens -
    A_DWN((ic-2)*4+2+7,(ic-1)*4+4,:)=mic*-u1(2,:,ic);                       %Amplitude longi sens -
    
    if ic<(ncapas-1)
        %U2 EN h(ic)
        A_DWN((ic-2)*4+2+8,(ic-1)*4+1,:)=mic* u2(1,:,ic).*exp1;                 %Amplitude trans sens +
        A_DWN((ic-2)*4+2+8,(ic-1)*4+2,:)=mic* u2(2,:,ic).*exp2;                 %Amplitude longi sens +
        A_DWN((ic-2)*4+2+8,(ic-1)*4+3,:)=mic*-u2(1,:,ic);                       %Amplitude trans sens -
        A_DWN((ic-2)*4+2+8,(ic-1)*4+4,:)=mic* u2(2,:,ic);                       %Amplitude longi sens -
    end
end

for ic=ncapas
    MAT=para.reg(1).sub(ic);
    
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    for i=1
        S11(i,:)= -1i*(k1(i,:,ic).*MAT.Ci(1,1).*u1(i,:,ic) +k2.*MAT.Ci(1,2).*u2(i,:,ic));
        S12(i,:)= -1i*MAT.Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
    end
    mic	= (-1)^ic;
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    %SIGMA11 EN 0
    A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)=mic* S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    %         A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)=mic* S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    
    %SIGMA12 EN 0
    A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    %         A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    
    %U1 EN 0
    A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);  	%Amplitude trans sens +
    %         A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    
    %         %U2 EN 0
    %         A_DWN((ic-2)*4+2+4,(ic-1)*4+1,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    %         A_DWN((ic-2)*4+2+4,(ic-1)*4+2,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul du vecteur source et resolution des coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT = para.reg(1).sub;

% le calcul du vecteur source se fait en utilisant
% la decomposition en ondes planes
%
% htop=0;
% for inch=1:ics-1
%     htop=htop+MAT(inch).h;
% end
%
% hbot=0;
% for inch=1:ics
%     hbot=hbot+MAT(inch).h;
% end

%%%%%%%%%%%%%%%%
% interface SE %
%%%%%%%%%%%%%%%%

%indice de los terminos en el vector fuente
is110   =(ics-2)*4+2+1;%SIGMA11 EN 0
is120   =(ics-2)*4+2+2;%SIGMA12 EN 0
iu10    =(ics-2)*4+2+3;%U1 EN 0
% iu20    =(ics-2)*4+2+4;%U2 EN 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues a l OP incidente %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%on suppose une onde transverse d amplitude unitaire de k1 imaginaire

%incidencia de ondas P
k1f     = k1(2,:,ics);%k1P
u1inc  	= u1(2,:,ics);
u2inc   = u2(2,:,ics);

u11     =-1i*k1f.*u1inc;
u21     =-1i*k1f.*u2inc;
u12     =-1i*k2 .*u1inc;
u22     =-1i*k2 .*u2inc;

Ci      = MAT(ics).Ci;
mic     =-(-1)^ics;

Fsource(is110,:)=Ci(1,1)*u11 + Ci(1,2)*u22;                                 %sigma11 en 0
Fsource(is120,:)=Ci(6,6)*(u12+u21);                                         %sigma12 en 0
if ics>1
    Fsource(iu10,:)=u1inc;                                                  %U1 EN 0
    %     Fsource(iu20,:)=u2;                                                	%U2 EN 0
end
Fsource=Fsource*mic;

%Resolution du systeme des conditions aux limites
% Solution   = A_DWN\Fsource;
Solution=zeros(4*(ncapas-1)+1,nk2);
for i=1:nk2
    Solution(:,i)   = A_DWN(:,:,i)\Fsource(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalisation en energie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%l origine de l'onde ds le SE doit etre au niveau de l interface
%SE/derniere couche pour que la normalisation soit correcte
E = 0;
% Ech=0;
% n=1000;
for icr=1:ncapas-1
    %E1
    
    h   = MAT(icr).h;
    kzS = k1(1,:,icr);
    kzP = k1(2,:,icr);
    
    %     %verification numerique de l integration
    %     dz=h/n;
    %     zrg=dz/2:dz:h;
    %     S1  = Solution(1+(icr-1)*4,:);
    %     S2  = Solution(2+(icr-1)*4,:);
    %     S3  = Solution(3+(icr-1)*4,:);
    %     S4  = Solution(4+(icr-1)*4,:);
    %
    %     for z0=zrg
    %         exp1= exp(-1i.*kzS.*z0);
    %         exp2= exp(-1i.*kzP.*z0);
    %         exp3= exp( 1i.*kzS.*(z0-h));
    %         exp4= exp( 1i.*kzP.*(z0-h));
    %
    %         U1KW = S1.*u1(1,:,icr).* exp1 + S2.*u1(2,:,icr).* exp2 + S3.*u1(1,:,icr).* exp3 - S4.*u1(2,:,icr).* exp4;
    %         U2KW = S1.*u2(1,:,icr).* exp1 + S2.*u2(2,:,icr).* exp2 - S3.*u2(1,:,icr).* exp3 + S4.*u2(2,:,icr).* exp4;
    %         Ech=Ech+sum(U1KW.*conj(U1KW)+U2KW.*conj(U2KW))*dz;
    %     end
    
    S1  = Solution(1+(icr-1)*4,:).*u1(1,:,icr);
    S2  = Solution(2+(icr-1)*4,:).*u1(2,:,icr);
    S3  = Solution(3+(icr-1)*4,:).*u1(1,:,icr);
    S4  =-Solution(4+(icr-1)*4,:).*u1(2,:,icr);
    
    E= E + 1i*(...
        + (S2*conj(S4) + S4*conj(S2))*kzP/(kzP^2 + abs(kzP)^2)*(exp(-h*kzP*1i) - exp(h*abs(kzP)^2*1i/kzP)) ...
        + (S1*conj(S3) + S3*conj(S1))*kzS/(kzS^2 + abs(kzS)^2)*(exp(-h*kzS*1i) - exp(h*abs(kzS)^2*1i/kzS)) ...
        + (S2*conj(S2) + S4*conj(S4))*kzP/(kzP^2 - abs(kzP)^2)*(exp(1i*h*(abs(kzP)^2/kzP-kzP)) - 1) ...
        + (S1*conj(S1) + S3*conj(S3))*kzS/(kzS^2 - abs(kzS)^2)*(exp(1i*h*(abs(kzS)^2/kzS-kzS)) - 1) ...
        - (S2*conj(S1) + S4*conj(S3))    /(kzP - conj(kzS))   *(1 - exp(1i*h*(conj(kzS)-kzP))) ...
        - (S1*conj(S2) + S3*conj(S4))    /(kzS - conj(kzP))   *(1 - exp(1i*h*(conj(kzP)-kzS))) ...
        - (S2*conj(S3) + S4*conj(S1))    /(kzP + conj(kzS))   *(exp(1i*h*conj(kzS)) - exp(-1i*h*kzP)) ...
        - (S1*conj(S4) + S3*conj(S2))    /(kzS + conj(kzP))   *(exp(1i*h*conj(kzP)) - exp(-1i*h*kzS)));
    
    S1  = Solution(1+(icr-1)*4,:).*u2(1,:,icr);
    S2  = Solution(2+(icr-1)*4,:).*u2(2,:,icr);
    S3  =-Solution(3+(icr-1)*4,:).*u2(1,:,icr);
    S4  = Solution(4+(icr-1)*4,:).*u2(2,:,icr);
    
    E= E + 1i*(...
        + (S2*conj(S4) + S4*conj(S2))*kzP/(kzP^2 + abs(kzP)^2)*(exp(-h*kzP*1i) - exp(h*abs(kzP)^2*1i/kzP)) ...
        + (S1*conj(S3) + S3*conj(S1))*kzS/(kzS^2 + abs(kzS)^2)*(exp(-h*kzS*1i) - exp(h*abs(kzS)^2*1i/kzS)) ...
        + (S2*conj(S2) + S4*conj(S4))*kzP/(kzP^2 - abs(kzP)^2)*(exp(1i*h*(abs(kzP)^2/kzP-kzP)) - 1) ...
        + (S1*conj(S1) + S3*conj(S3))*kzS/(kzS^2 - abs(kzS)^2)*(exp(1i*h*(abs(kzS)^2/kzS-kzS)) - 1) ...
        - (S2*conj(S1) + S4*conj(S3))    /(kzP - conj(kzS))   *(1 - exp(1i*h*(conj(kzS)-kzP))) ...
        - (S1*conj(S2) + S3*conj(S4))    /(kzS - conj(kzP))   *(1 - exp(1i*h*(conj(kzP)-kzS))) ...
        - (S2*conj(S3) + S4*conj(S1))    /(kzP + conj(kzS))   *(exp(1i*h*conj(kzS)) - exp(-1i*h*kzP)) ...
        - (S1*conj(S4) + S3*conj(S2))    /(kzS + conj(kzP))   *(exp(1i*h*conj(kzP)) - exp(-1i*h*kzS)));
end
icr = ncapas;

% n=10000;
% h=10;
% dz=h/n;
% zrg=dz/2:dz:h;
% S1  = 1                      ;
% S2  = Solution(1+(icr-1)*4,:);
% for z0=zrg
%     exp1= exp(-1i.*k1(1,:,icr).*z0);
%     exp2= exp(-1i.*k1(2,:,icr).*z0);
%
%     U1KW = S1.*u1(1,:,icr).* exp1 + S2.*u1(2,:,icr).* exp2;
%     U2KW = S1.*u2(1,:,icr).* exp1 + S2.*u2(2,:,icr).* exp2;
%     Ech=Ech+sum(U1KW.*conj(U1KW)+U2KW.*conj(U2KW))*dz;
% end%OK l integration numerique donne la meme chose que l analytique

kzS =-imag(k1(1,:,icr));
kzP =-imag(k1(2,:,icr));
% U1KW(1:nk2) = 1.*u1(1,:,icr).* exp1 + Solution(1+(icr-1)*4,:).*u1(2,:,icr).* exp2;
% U2KW(1:nk2) = 1.*u2(1,:,icr).* exp1 + Solution(1+(icr-1)*4,:).*u2(2,:,icr).* exp2;
S1  = Solution(1+(icr-1)*4,:) .*u1(1,:,icr);
S2  = 1                       .*u1(2,:,icr);

E   = E + abs(S2)^2/(2*kzP) + abs(S1)^2/(2*kzS) + (S1*conj(S2)+S2*conj(S1))/(kzP + kzS);

S1  = Solution(1+(icr-1)*4,:) .*u2(1,:,icr);
S2  = 1                       .*u2(2,:,icr);

E   = E + abs(S2)^2/(2*kzP) + abs(S1)^2/(2*kzS) + (S1*conj(S2)+S2*conj(S1))/(kzP + kzS);

E   = real(sqrt(E));%pour reinjecter ds les deplacements

%normalisation des amplitudes
Solution=Solution/E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des champs au niveau du recepteur %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for irz=1:nrecz
    ixr     = find(izr0==irz);
    nrecx   = length(ixr);
    
    goU     = sum(salu(ixr(1:nrecx)));
    goS     = sum(sals(ixr(1:nrecx)));
    %Ecriture de la solution en fonction de la position du recepteur
    zricr   = zr0(irz);
    
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
            U1KW = Solution(1+(icr-1)*4,:).*u1(2,:,icr).* exp2;
            U2KW = Solution(1+(icr-1)*4,:).*u2(2,:,icr).* exp2;
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
            S11KW = Solution(1+(icr-1)*4,:).*S11(2,:).* exp2;
            S22KW = Solution(1+(icr-1)*4,:).*S22(2,:).* exp2;
            S12KW = Solution(1+(icr-1)*4,:).*S12(2,:).* exp2;
        end
        S11KW(isnan(S11KW))=0;
        S22KW(isnan(S22KW))=0;
        S12KW(isnan(S12KW))=0;
    end
    
    %rajout des champs incidents
    if icr==ics
        zrs=zr0(irz)-zs;
        
        %incidencia de ondas S
        expt    = exp(-1i.*k1f.*zrs);
        u1inc 	= u1(1,:,ics)./E.*expt;
        u2inc 	= u2(1,:,ics)./E.*expt;
        
        u11     =-1i*k1f.*u1inc;
        u21     =-1i*k1f.*u2inc;
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