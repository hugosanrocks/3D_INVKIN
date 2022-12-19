function DWN=calcul_A_DWN_PSV_Ncapas_HS(para,DWN)
% METHODE : @ da los valores de k1 (componente vertical de P y S) por un medio isotropo
%           @ da los vectores propio, polarizacion associados a k1 y k2
%           @ calculo de la matrice de los coefficiente de la sol homogenea
%             a partir de las condiciones en las interfases

ncapas      = para.nsubmed;
k2          = DWN.k2;
nk2       	= length(k2(1,:));
if nk2==1
    if isfield(DWN,'omegac')
        nk2=length(DWN.omegac);
    end
end
u1          = zeros(2,nk2,ncapas);
u2          = zeros(2,nk2,ncapas);
k1          = zeros(2,nk2,ncapas);

S11         = zeros(2,nk2);
S12         = zeros(2,nk2);

A_DWN       = zeros(4*(para.nsubmed-1)+2,4*(para.nsubmed-1)+2,nk2);
% detA        = zeros(nk2,1);

%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
for ic=1:ncapas
    ksi         = para.reg(1).sub(ic).ksi;
    kpi         = para.reg(1).sub(ic).kpi;
    k1(1,:,ic)  = sqrt(ksi.^2-k2.^2); %k1S
    k1(2,:,ic)  = sqrt(kpi.^2-k2.^2); %k1P
end
k1 =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));
for ic=1:ncapas
    u1(1,:,ic)  =  k2;          % S
    u2(1,:,ic)  = -k1(1,:,ic);  % S
    u1(2,:,ic)  =  k1(2,:,ic);  % P
    u2(2,:,ic)  =  k2;          % P
end
un      = sqrt(abs(u1).^2+abs(u2).^2);
u1      = u1./un;
u2      = u2./un;
DWN.k1  = k1;
DWN.u1  = u1;
DWN.u2  = u2;

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

    if ncapas>1
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
    
    %U2 EN h(ic)
    A_DWN((ic-2)*4+2+8,(ic-1)*4+1,:)=mic* u2(1,:,ic).*exp1;                 %Amplitude trans sens +
    A_DWN((ic-2)*4+2+8,(ic-1)*4+2,:)=mic* u2(2,:,ic).*exp2;                 %Amplitude longi sens +
    A_DWN((ic-2)*4+2+8,(ic-1)*4+3,:)=mic*-u2(1,:,ic);                       %Amplitude trans sens -
    A_DWN((ic-2)*4+2+8,(ic-1)*4+4,:)=mic* u2(2,:,ic);                       %Amplitude longi sens -
end

if ncapas>1
    for ic=ncapas
        MAT=para.reg(1).sub(ic);
        
        %***********************************************************%
        %   calcul de la contrainte issue de la solution homogene   %
        %***********************************************************%
        for i=1:2
            S11(i,:)= -1i*(k1(i,:,ic).*MAT.Ci(1,1).*u1(i,:,ic) +k2.*MAT.Ci(1,2).*u2(i,:,ic));
            S12(i,:)= -1i*MAT.Ci(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
        end
        mic	= (-1)^ic;
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        %SIGMA11 EN 0
        A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)=mic* S11(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
        A_DWN((ic-2)*4+2+1,(ic-1)*4+2,:)=mic* S11(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
        
        %SIGMA12 EN 0
        A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* S12(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude trans sens +
        A_DWN((ic-2)*4+2+2,(ic-1)*4+2,:)=mic* S12(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude longi sens +
        
        %U1 EN 0
        A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);  	%Amplitude trans sens +
        A_DWN((ic-2)*4+2+3,(ic-1)*4+2,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
        
        %U2 EN 0
        A_DWN((ic-2)*4+2+4,(ic-1)*4+1,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude trans sens +
        A_DWN((ic-2)*4+2+4,(ic-1)*4+2,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude longi sens +
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % inversion de la matrice %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:nk2
%     A_DWN(:,:,i)=inv(A_DWN(:,:,i));
% end

DWN.A_DWN = A_DWN;