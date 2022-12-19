function [UXW,SXW]=calcul_US_DWN_PSV_Ncapas_HS_incOPI(para,xr,zr0,izr0,salu,sals,kxk)
% caso de incidencia de ondas inhomogenea en multi estratos
% en este caso se supone por lo menos 1 estrato sobre un semi espacio, el
% caso del semi-espacio esta procesado a parte pero se puede poner estratos
% de misma propiedades para verificar
% Se impone una incidencia de amplitud unitaria
% la onda incidente se impone en el estrato de uno de los receptores para
% evitar problemas numericos
% se busca las amplitudes de las ondas en los estratos tomando en cuenta la
% incidencia de la onda que se impuso
% asi el sistema se reduce de una dimension con respeto al caso normal
% se quita la linea s11 donde esta la fuente, y la amplitud de la onda inc
% se normaliza la energia de cada modo de manera que la integral con
% respecto a la profundidad es unitaria


ncapas      = para.nsubmed;
nk2         = length(kxk);
nrec        = length(xr);
nrecz       = length(zr0);

%identificacion del estrato del receptor de minima profundidad
zrm         = min(zr0);
ics         = 1;
MAT         = para.reg(1).sub;
hc          = MAT(1).h;
while zrm>hc && ics<ncapas
    ics     = ics+1;
    hc      = hc+MAT(ics).h;
end
% ics=ncapas+1;

% se considera que la onda incidente proviene del estrato del receptor de
% minima profundidad
xs          = mean(xr);
nxs       	= length(xs); %### a desarollar para varias incidencias

%inicialisaciones
UXW         = zeros(2  ,nrec,nxs);
SXW         = zeros(2,2,nrec,nxs);
S11         = zeros(2,nk2);
S12        	= zeros(2,nk2);
S22         = zeros(2,nk2);

u1          = zeros(2,nk2,ncapas);
u2      	= zeros(2,nk2,ncapas);

%se reduce una dimension c/r al caso normal
%se escogio suprimir u2 en el ultimo estrato (la ultima ecuacion)
%y la amplitud de la onda que se toma como incidente (S bajando)
Fsource     = zeros(4*(ncapas-1)+1,nk2);
Solution     =zeros(4*(ncapas-1)+1,nk2);
A_DWN     	= zeros(4*(para.nsubmed-1)+1,4*(para.nsubmed-1)+1,nk2);


% identification des composantes du nombre d onde incident
% l incidence provient du semi-espace ### quizas a cambiar si incidencia en
% otros estratos
ksi         = para.reg(1).sub(ncapas).ksi;
k2          = kxk*ksi;

%**********************************************************************%
%  calcul des nombres d'onde et polarisations de la solution homogene  %
%**********************************************************************%
k1              = zeros(2,nk2,ncapas);
for ic=1:ncapas
    ksi         = para.reg(1).sub(ic).ksi;
    kpi         = para.reg(1).sub(ic).kpi;
    k1(1,:,ic)  = sqrt(ksi^2-k2.^2); %k1S
    k1(2,:,ic)  = sqrt(kpi^2-k2.^2); %k1P
end
k1              =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));

for ic=1:ncapas
    u1(1,:,ic)  =  k2;          % S
    u2(1,:,ic)  = -k1(1,:,ic);  % S
    u1(2,:,ic)  =  k1(2,:,ic);  % P
    u2(2,:,ic)  =  k2;          % P
end

%normalisacion de los desplazamientos
un              = sqrt(abs(u1).^2+abs(u2).^2);
u1              = u1./un;
u2              = u2./un;

%************************************************%
% calcul de la matriz de condicion de interfaces %
%************************************************%

%les ondes en exp(-1i.*k1(:,ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*k1(:,ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    h   = MAT(ic).h;
    Ci  = MAT(ic).Ci;
    mic	= (-1)^ic;
    %translation d indice pour assurer la construction de la matrice sans
    %prendre en compte l incidence imposee
    itr = real(ics<=ic);
    
    exp1= exp(-1i.*k1(1,:,ic).*h);      %condition en h
    exp2= exp(-1i.*k1(2,:,ic).*h);
    exp3= exp( 1i.*k1(1,:,ic).*(0-h));  %condition en 0
    exp4= exp( 1i.*k1(2,:,ic).*(0-h));
    
    for i=1:2
        S11(i,:)= -1i*(k1(i,:,ic).*Ci(1,1).*u1(i,:,ic) +k2.*Ci(1,2).*u2(i,:,ic));   %sigma_11 %Amplitude longi & transverse sens z + (profondeur)
        S12(i,:)= -1i*Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));            %sigma_12 %Amplitude longi & transverse sens z +
    end
    % S11(1,-)=-S11(1,+) pour k1-=-k1+ & S
    % S12(1,-)= S12(1,+)
    % S11(2,-)= S11(2,+) pour k1-=-k1+ & P
    % S12(2,-)=-S12(2,+)
    
    %SIGMA11 EN 0
    if ics~=ic
        A_DWN(1,1,:)= mic*S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);                 %Amplitude trans sens +
        A_DWN(1,2,:)= mic*S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);                	%Amplitude longi sens +
        A_DWN(1,3,:)=-mic*S11(1,:).*exp3;                                       %Amplitude trans sens -
        A_DWN(1,4,:)= mic*S11(2,:).*exp4;                                       %Amplitude longi sens -
    end
    
    %SIGMA12 EN 0
    if ics~=ic
        A_DWN(2,1,:)= mic*S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);                 %Amplitude trans sens +
    end
    A_DWN(2-itr,2-itr,:)= mic*S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);                 %Amplitude longi sens +
    A_DWN(2-itr,3-itr,:)= mic*S12(1,:).*exp3;                                       %Amplitude trans sens -
    A_DWN(2-itr,4-itr,:)=-mic*S12(2,:).*exp4;                                       %Amplitude longi sens -
    
    %SIGMA11 EN h(ic)
    if ics~=ic+1
        if ics~=ic
            A_DWN(3-itr,1,:)= mic*S11(1,:).*exp1;                                       %Amplitude trans sens +
        end
        A_DWN(3-itr,2-itr,:)= mic*S11(2,:).*exp2;                                       %Amplitude longi sens +
        A_DWN(3-itr,3-itr,:)=-mic*S11(1,:);                                             %Amplitude trans sens -
        A_DWN(3-itr,4-itr,:)= mic*S11(2,:);                                             %Amplitude longi sens -
        itr1                = itr;
    else
        itr1                = itr+1;
    end
    
    %SIGMA12 EN h(ic)
    if ics~=ic
        A_DWN(4-itr1,1,:) =mic*S12(1,:).*exp1;                                       %Amplitude trans sens +
    end
    A_DWN(4-itr1,2-itr,:)= mic*S12(2,:).*exp2;                                       %Amplitude longi sens +
    A_DWN(4-itr1,3-itr,:)= mic*S12(1,:);                                             %Amplitude trans sens -
    A_DWN(4-itr1,4-itr,:)=-mic*S12(2,:);                                             %Amplitude longi sens -
    
    %U1 EN h(ic)
    if ics~=ic
        A_DWN(5-itr1,1,:)= mic*u1(1,:,ic).*exp1;                                     %Amplitude trans sens +
    end
    A_DWN(5-itr1,2-itr,:)= mic*u1(2,:,ic).*exp2;                                     %Amplitude longi sens +
    A_DWN(5-itr1,3-itr,:)= mic*u1(1,:,ic);                                           %Amplitude trans sens -
    A_DWN(5-itr1,4-itr,:)=-mic*u1(2,:,ic);                                           %Amplitude longi sens -
    
    %U2 EN h(ic)
    if ics~=ic
        A_DWN(6-itr1,1,:)= mic*u2(1,:,ic).*exp1;                                 %Amplitude trans sens +
    end
    A_DWN(6-itr1,2-itr,:)= mic*u2(2,:,ic).*exp2;                                 %Amplitude longi sens +
    A_DWN(6-itr1,3-itr,:)=-mic*u2(1,:,ic);                                       %Amplitude trans sens -
    A_DWN(6-itr1,4-itr,:)= mic*u2(2,:,ic);                                       %Amplitude longi sens -
end

for ic=2:(ncapas-1)
    h   = MAT(ic).h;
    Ci  = MAT(ic).Ci;
    itr = real(ics<=ic);
    
    for i=1:2
        S11(i,:)= -1i*(k1(i,:,ic).*Ci(1,1).*u1(i,:,ic) +k2.*Ci(1,2).*u2(i,:,ic));
        S12(i,:)= -1i*Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
    end
    
    exp1= exp(-1i.*k1(1,:,ic).*h);
    exp2= exp(-1i.*k1(2,:,ic).*h);
    exp3= exp( 1i.*k1(1,:,ic).*(0-h));
    exp4= exp( 1i.*k1(2,:,ic).*(0-h));
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+1-itr,:)=mic* S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+2-itr,:)=mic* S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+3-itr,:)=mic*-S11(1,:).*exp3;                           %Amplitude trans sens -
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+4-itr,:)=mic* S11(2,:).*exp4;                           %Amplitude longi sens -
    end
    %SIGMA12 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+1-itr,:)=mic* S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+2-itr,:)=mic* S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+3-itr,:)=mic* S12(1,:).*exp3;                           %Amplitude trans sens -
    A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+4-itr,:)=mic*-S12(2,:).*exp4;                           %Amplitude longi sens -
    
    %U1 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+1-itr,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+2-itr,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+3-itr,:)=mic* u1(1,:,ic).*exp3;                       	%Amplitude trans sens -
    A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+4-itr,:)=mic*-u1(2,:,ic).*exp4;                         %Amplitude longi sens -
    
    %U2 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+1-itr,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+2-itr,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+3-itr,:)=mic*-u2(1,:,ic).*exp3;                         %Amplitude trans sens -
    A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+4-itr,:)=mic* u2(2,:,ic).*exp4;                         %Amplitude longi sens -
    
    %SIGMA11 EN h(ic)
    if ics~=ic+1
        if ics~=ic
            A_DWN((ic-2)*4+2+5-itr,(ic-1)*4+1-itr,:)=mic* S11(1,:).*exp1;                   %Amplitude trans sens +
        end
        A_DWN((ic-2)*4+2+5-itr,(ic-1)*4+2-itr,:)=mic* S11(2,:).*exp2;                   %Amplitude longi sens +
        A_DWN((ic-2)*4+2+5-itr,(ic-1)*4+3-itr,:)=mic*-S11(1,:);                         %Amplitude trans sens -
        A_DWN((ic-2)*4+2+5-itr,(ic-1)*4+4-itr,:)=mic* S11(2,:);                         %Amplitude longi sens -
        itr1                = itr;
    else
        itr1                = itr+1;
    end
    
    %SIGMA12 EN h(ic)
    if ics~=ic
        A_DWN((ic-2)*4+2+6-itr1,(ic-1)*4+1-itr,:)=mic* S12(1,:).*exp1;                   %Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+6-itr1,(ic-1)*4+2-itr,:)=mic* S12(2,:).*exp2;                   %Amplitude longi sens +
    A_DWN((ic-2)*4+2+6-itr1,(ic-1)*4+3-itr,:)=mic* S12(1,:);                         %Amplitude trans sens -
    A_DWN((ic-2)*4+2+6-itr1,(ic-1)*4+4-itr,:)=mic*-S12(2,:);                         %Amplitude longi sens -
    
    %U1 EN h(ic)
    if ics~=ic
        A_DWN((ic-2)*4+2+7-itr1,(ic-1)*4+1-itr,:)=mic* u1(1,:,ic).*exp1;                 %Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+7-itr1,(ic-1)*4+2-itr,:)=mic* u1(2,:,ic).*exp2;                 %Amplitude longi sens +
    A_DWN((ic-2)*4+2+7-itr1,(ic-1)*4+3-itr,:)=mic* u1(1,:,ic);                       %Amplitude trans sens -
    A_DWN((ic-2)*4+2+7-itr1,(ic-1)*4+4-itr,:)=mic*-u1(2,:,ic);                       %Amplitude longi sens -
    
    %U2 EN h(ic)
    if ics~=ic
        A_DWN((ic-2)*4+2+8-itr1,(ic-1)*4+1-itr,:)=mic* u2(1,:,ic).*exp1;                 %Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+8-itr1,(ic-1)*4+2-itr,:)=mic* u2(2,:,ic).*exp2;                 %Amplitude longi sens +
    A_DWN((ic-2)*4+2+8-itr1,(ic-1)*4+3-itr,:)=mic*-u2(1,:,ic);                       %Amplitude trans sens -
    A_DWN((ic-2)*4+2+8-itr1,(ic-1)*4+4-itr,:)=mic* u2(2,:,ic);                       %Amplitude longi sens -
end

for ic=ncapas
    Ci  = MAT(ic).Ci;
    itr = real(ics<=ic);
    
    %   calcul de la contrainte issue de la solution homogene   %
    for i=2
        S11(i,:)= -1i*(k1(i,:,ic).*Ci(1,1).*u1(i,:,ic) +k2.*Ci(1,2).*u2(i,:,ic));
        S12(i,:)= -1i*Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
    end
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+1-itr,:)=mic* S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
        A_DWN((ic-2)*4+2+1-itr,(ic-1)*4+2-itr,:)=mic* S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    end
    
    %SIGMA12 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+1-itr,:)=mic* S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+2-itr,(ic-1)*4+2-itr,:)=mic* S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
    
    %U1 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+1-itr,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);  	%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+3-itr,(ic-1)*4+2-itr,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    
    %U2 EN 0
    if ics~=ic
        A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+1-itr,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
    end
    A_DWN((ic-2)*4+2+4-itr,(ic-1)*4+2-itr,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
end

% if ics==ncapas+1
%     %solution promedia
%     %el determinante es cercano a cero pero no lo es
%     Am1=inv(A_DWN);
%     %normalizacion del primer termino
%     NA=(4*(ncapas-1)+2);
%     for i=1:NA
%         Am1(:,i)=Am1(:,i)/Am1(1,i);
%     end
%     [~,ind]=min(sum(abs(Am1).^2,1));
%     Solution=Am1(:,4);%vecteur propre
% %     Solution=mean(Am1,2);%vecteur propre
% for kkk=1:NA
%     sum(abs(A_DWN*Am1(:,i)).^2)
% end
% else
    
    %**************************%
    % calcul du vecteur source %
    %**************************%
    
    %indice de los terminos en el vector fuente
    %se suprima la ultima linea en caso de que la incidencia viene del semi espacio
    if ics==1
        %is110=1;           %SIGMA11 EN 0
        is120=1;            %SIGMA12 EN 0
        is11h=2;            %SIGMA11 EN h
        is12h=3;            %SIGMA12 EN h
        iu1h =4;            %U1 EN h
        iu2h =5;            %U2 EN h
    elseif ics==ncapas
        %is110=(ics-2)*4+2+1;   %SIGMA11 EN 0
        is120=(ics-2)*4+2+1;    %SIGMA12 EN 0
        iu10 =(ics-2)*4+2+2;    %U1 EN 0
        iu20 =(ics-2)*4+2+3;    %U2 EN 0
    else %1<ics<ncapas
        %is110=(ics-2)*4+2+1;   %SIGMA11 EN 0
        is120=(ics-2)*4+2+1;    %SIGMA12 EN 0
        iu10 =(ics-2)*4+2+2;    %U1 EN 0
        iu20 =(ics-2)*4+2+3;    %U2 EN 0
        is11h=(ics-2)*4+2+4;    %SIGMA11 EN h
        is12h=(ics-2)*4+2+5;    %SIGMA12 EN h
        iu1h =(ics-2)*4+2+6;    %U1 EN h
        iu2h =(ics-2)*4+2+7;    %U2 EN h
    end
    
    % htop=0;
    % for inch=1:ics-1
    %     htop=htop+MAT(inch).h;
    % end
    %
    % hbot=0;
    % for inch=1:ics
    %     hbot=hbot+MAT(inch).h;
    % end
    
    
    % calcul des conditions aux interface dues a l OP incidente %
    % on suppose une onde transverse d amplitude unitaire de k1 imaginaire en
    % direction vers les z + (vers la profondeur)
    
    %incidencia de ondas S
    k1f     = k1(1,:,ics);%k1S
    u1inc  	= u1(1,:,ics);
    u2inc   = u2(1,:,ics);
    
    u11     =-1i*k1f.*u1inc;
    u21     =-1i*k1f.*u2inc;
    u12     =-1i*k2 .*u1inc;
    u22     =-1i*k2 .*u2inc;
    Ci      = MAT(ics).Ci;
    S11     = Ci(1,1)*u11 + Ci(1,2)*u22;
    S12     = Ci(6,6)*(u12+u21);
    if ics<ncapas
        h       = MAT(ics).h;
        exp1    = exp(-1i.*k1f*h);
    end
    
    
    % calcul du terme source %
    % Fsource(is110,:)=S11;                       %sigma11 en 0
    Fsource(is120,:)=S12;                       %sigma12 en 0
    if ics>1
        Fsource(iu10,:)=u1inc;                  %u1 EN 0
        Fsource(iu20,:)=u2inc;              %u2 EN 0
    end
    if ics<ncapas
        Fsource(is11h,:)=S11  .*exp1;       	%sigma11 en h
        Fsource(is12h,:)=S12  .*exp1;         	%sigma12 en h
        Fsource(iu1h ,:)=u1inc.*exp1;           %U1 EN h
        Fsource(iu2h ,:)=u2inc.*exp1;           %U2 EN h
    end
    mic     =-(-1)^ics;
    Fsource = Fsource*mic;
    
    %*****************************%
    % resolution des coefficients %
    %*****************************%
    
    %Resolution du systeme des conditions aux limites
    % Solution   = A_DWN\Fsource;
    for i=1:nk2
        Solution(:,i)   = A_DWN(:,:,i)\Fsource(:,i);
    end
    
    %on reinsere la composante incidente ds la solution
    %et on supprime plus tard la contribution du champ incident
    n0 = 4*(para.nsubmed-1)+1;
    for i=1:nk2
        Solution((((ics-1)*4+1):n0)+1,i) = Solution(((ics-1)*4+1):n0,i);
    end
    Solution((ics-1)*4+1,i)   = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalisation en energie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalisation dans les couches
E = 0;
for ic=1:ncapas-1
    h   = MAT(ic).h;
    kzS = k1(1,:,ic);
    kzP = k1(2,:,ic);
    
    %E(u1)
    S1  = Solution(1+(ic-1)*4,:).*u1(1,:,ic);%Amplitude trans sens +
    S2  = Solution(2+(ic-1)*4,:).*u1(2,:,ic);%Amplitude longi sens +
    S3  = Solution(3+(ic-1)*4,:).*u1(1,:,ic);%Amplitude trans sens -
    S4  =-Solution(4+(ic-1)*4,:).*u1(2,:,ic);%Amplitude longi sens -
    
    E   = E + NormE_PSV(kzS,kzP,S1,S2,S3,S4,h);
    
    %E(u2)
    S1  = Solution(1+(ic-1)*4,:).*u2(1,:,ic);
    S2  = Solution(2+(ic-1)*4,:).*u2(2,:,ic);
    S3  =-Solution(3+(ic-1)*4,:).*u2(1,:,ic);
    S4  = Solution(4+(ic-1)*4,:).*u2(2,:,ic);
    
    E   = E + NormE_PSV(kzS,kzP,S1,S2,S3,S4,h);
end


% ic=1;
%     h   = MAT(ic).h;
%     kzS = k1(1,:,ic);
%     kzP = k1(2,:,ic);
% %verification numerique de l integration
% n   = 1000;
% dz  = h/n;
% zrg = dz/2:dz:h;
% S1  = Solution(1+(ic-1)*4,:);
% S2  = Solution(2+(ic-1)*4,:);
% S3  = Solution(3+(ic-1)*4,:);
% S4  = Solution(4+(ic-1)*4,:);
%
% exp1= exp(-1i.*kzS.*zrg);
% exp2= exp(-1i.*kzP.*zrg);
% exp3= exp( 1i.*kzS.*(zrg-h));
% exp4= exp( 1i.*kzP.*(zrg-h));
%
% U1KW = S1.*u1(1,:,ic).* exp1 + S2.*u1(2,:,ic).* exp2 + S3.*u1(1,:,ic).* exp3 - S4.*u1(2,:,ic).* exp4;
% U2KW = S1.*u2(1,:,ic).* exp1 + S2.*u2(2,:,ic).* exp2 - S3.*u2(1,:,ic).* exp3 + S4.*u2(2,:,ic).* exp4;
% E=sum(U1KW.*conj(U1KW)+U2KW.*conj(U2KW))*dz;



%normalisation dans le SE
ic  = ncapas;

% n=10000;
% h=10;
% dz=h/n;
% zrg=dz/2:dz:h;
% S1  = Solution(1+(ic-1)*4,:);
% S2  = Solution(2+(ic-1)*4,:);
% exp1= exp(-1i.*k1(1,:,ic).*zrg);
% exp2= exp(-1i.*k1(2,:,ic).*zrg);
%
% U1KW = S1.*u1(1,:,ic).* exp1 + S2.*u1(2,:,ic).* exp2;
% U2KW = S1.*u2(1,:,ic).* exp1 + S2.*u2(2,:,ic).* exp2;
% Ech  = sum(U1KW.*conj(U1KW)+U2KW.*conj(U2KW))*dz;
% %OK l integration numerique donne la meme chose que l analytique

kzS     =-imag(k1(1,:,ic));
kzP     =-imag(k1(2,:,ic));

S1      = Solution(1+(ic-1)*4,:).*u1(1,:,ic);
S2      = Solution(2+(ic-1)*4,:).*u1(2,:,ic);
E       = E + abs(S2)^2/(2*kzP) + abs(S1)^2/(2*kzS) + (S1*conj(S2)+S2*conj(S1))/(kzP + kzS);

S1      = Solution(1+(ic-1)*4,:).*u2(1,:,ic);
S2      = Solution(2+(ic-1)*4,:).*u2(2,:,ic);
E       = E + abs(S2)^2/(2*kzP) + abs(S1)^2/(2*kzS) + (S1*conj(S2)+S2*conj(S1))/(kzP + kzS);

E       = real(sqrt(E));%pour reinjecter ds les deplacements

%normalisation des amplitudes
Solution= Solution/E;

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
            U1KW = ...
                +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2;
            U2KW = ...
                +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2;
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