function DWN=calcul_A_DWN_SH_Ncapas_HS(para,DWN)
% METHODE : @ da los valores de k1 (componente vertical de P y S) por un medio isotropo
%           @ da los vectores propio, polarizacion associados a k1 y k2
%           @ calculo de la matrice de los coefficiente de la sol homogenea
%             a partir de las condiciones en las interfases

% ! ojo aqui el eje 2 es el eje horizontal y 1 el vertical, es diferente a
% las convenciones en el IBEM
ncapas      = para.nsubmed;
k2          = DWN.k2;
nk2       	= length(k2(1,:));
if nk2==1
    nk2=length(DWN.omegac);
end
k1          = zeros(nk2,ncapas);
A_DWN       = zeros(2*(para.nsubmed-1)+1,2*(para.nsubmed-1)+1,nk2);

%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
for ic=1:ncapas
    ksi 	= para.reg(1).sub(ic).ksi;
    k1(:,ic)= sqrt(ksi.^2-k2.^2); %k1S
end
k1 =-k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));

DWN.k1=k1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%les ondes en exp(-1i.*k1(:,ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*k1(:,ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    MAT = para.reg(1).sub(ic);
    S12 =-1i*k1(:,ic)*MAT.Ci(6,6); %sigma_zy sigma12
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    %SIGMA12 EN 0
    A_DWN(1,1,:)= mic*S12; %.*exp(-1i.*k1(:,ic).*0);                        %Amplitude trans sens +
    if ncapas>1
        exp2=exp(-1i.*k1(:,ic).*MAT.h);
        exp4=exp( 1i.*k1(:,ic).*(0-MAT.h));

        A_DWN(1,2,:)=-mic*S12.*exp4;                                        %Amplitude trans sens -
        
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
    S12 =-1i*k1(:,ic)*MAT.Ci(6,6); %sigma_12

    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    exp2= exp(-1i.*k1(:,ic).*MAT.h);
    exp4= exp( 1i.*k1(:,ic).*(0-MAT.h));
    mic	= (-1)^ic;
    
    %SIGMA12 EN 0
    A_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*S12; %.*exp(-1i.*k1(:,ic).*0);	%Amplitude trans sens +
    A_DWN((ic-2)*2+1+1,(ic-1)*2+2,:)=-mic*S12.*exp4; %.*exp( 1i.*k1(:,ic).*0);%Amplitude trans sens -
    
    %U2 EN 0
    A_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic; %.*exp(-1i.*k1(:,ic).*0);        %Amplitude trans sens +
    A_DWN((ic-2)*2+1+2,(ic-1)*2+2,:)= mic.*exp4; %.*exp( 1i.*k1(:,ic).*0);	%Amplitude trans sens -
    
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
        S12 =-1i*k1(:,ic)*MAT.Ci(6,6); %sigma_12
        
        mic	= (-1)^ic;
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        %SIGMA12 EN 0
        A_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*S12; %.*exp(-1i.*k1(:,ic).*0);		%Amplitude trans sens +
        %         A_DWN((ic-2)*2+1+1,(ic-1)*2+2,:)=mic* S12; %.*exp(-1i.*k1(:,ic).*0);	    %Amplitude trans sens -
        
        %U2 EN 0
        A_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic; %.*exp(-1i.*k1(:,ic).*0);	%Amplitude trans sens +
        %         A_DWN((ic-2)*2+1+2,(ic-1)*2+2,:)=mic; %.*exp(-1i.*k1(:,ic).*0);	%Amplitude trans sens -
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion de la matrice %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nk2
%     A_DWN(:,:,i)=inv(A_DWN(:,:,i));
% end
DWN.A_DWN = A_DWN;

% %verifier que P est tjs le meme
% [L,U,P,Q] = lu(A_DWN(:,:,5),'vector');
% s =  det(L);        
% det(A) = s*prod(diag(U));