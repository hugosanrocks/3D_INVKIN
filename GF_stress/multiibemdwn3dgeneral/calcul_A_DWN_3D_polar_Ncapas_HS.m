function DWN=calcul_A_DWN_3D_polar_Ncapas_HS(para,DWN)
% Matrices A,B del medio estratificado (DWN) solución en polares (rápido)
% METHODE : @ da los valores de kz (componente vertical de P,SV,SH) por un medio isotropo
%           @ da los vectores propio, polarizacion associados a kr y kz
%           @ calculo de la matrice de los coefficiente de la sol homogenea
%             a partir de las condiciones en las interfases
% OJO!!! se cambia los ejes con respectos a los modelos PSV 1=x, 2==y, 3==z
% tic
ncapas      = para.nsubmed;
kr          = DWN.kr; % (k2)
nkr       	= length(kr(1,:));
if nkr==1
    if isfield(DWN,'omegac')
        nkr=length(DWN.omegac);
    end
end
MAT         = para.reg(1).sub;

u1          = zeros(2,nkr,ncapas);
u2          = zeros(2,nkr,ncapas);
kz          = zeros(2,nkr,ncapas);
%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
for ic=1:ncapas
    ksi         = para.reg(1).sub(ic).ksi;
    kpi         = para.reg(1).sub(ic).kpi;
    kz(1,:,ic)  = sqrt(ksi.^2-kr.^2); %kzS nu
    kz(2,:,ic) 	= sqrt(kpi.^2-kr.^2); %kzP gamma
end
kz  =-kz.*(sign(imag(kz)).*(imag(kz)~=0)+(imag(kz)==0));
DWN.kz=kz; % (.k1)
for ic=1:ncapas
    u1(1,:,ic)  =  kr;          % S
    u2(1,:,ic)  = -kz(1,:,ic);  % S
    u1(2,:,ic)  =  kz(2,:,ic);  % P
    u2(2,:,ic)  =  kr;          % P
end
un      = sqrt(abs(u1).^2+abs(u2).^2);
u1      = u1./un;
u2      = u2./un;
DWN.u1  = u1;
DWN.u2  = u2;


% spmd
% codist = codistributor('1d',2); %codistrituir la variable sobre la dimensión 2
codist=1;
Szz         = zeros(3,nkr,codist);
Szr         = zeros(3,nkr,codist);
% Szt         = zeros(3,nkr);
Uz          = zeros(3,nkr,codist);
Ur          = zeros(3,nkr,codist);
% Ut          = zeros(3,nkr);
% codist = codistributor('1d',3); %codistrituir la variable sobre la dimensión 2
A_DWN       = zeros(4*(ncapas-1)+2,4*(ncapas-1)+2,nkr,codist);%P-SV
B_DWN       = zeros(2*(ncapas-1)+1,2*(ncapas-1)+1,nkr,codist);%SH
%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%les ondes en exp(-1i.*kz(:,ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*kz(:,ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    h  = para.reg(1).sub(ic).h;
    xi = kz(1,:,ic).^2-kr.^2;
    
    Szz(1,:)=-xi;                                                           %Amplitude Phi sens z +
    Szz(2,:)=-2i*kz(1,:,ic).*kr.^2;                                         %Amplitude Psi z +
    Szz(3,:)= 0;                                                            %Amplitude Xsi z +
    
    Szr(1,:)= 2i*kz(2,:,ic);                                                %Amplitude Phi z +
    Szr(2,:)= xi;                                                           %Amplitude Psi z +
    Szr(3,:)= 1i*kz(1,:,ic);                                                %Amplitude Xsi z +
    
    Szz     = MAT(ic).Ci(6,6)*Szz;
    Szr     = MAT(ic).Ci(6,6)*Szr;
    
    %     Szt(1,:)= 2i*kz(2,:,ic);                                                %Amplitude Phi z +
    %     Szt(2,:)= xi;                                                           %Amplitude Psi z +
    %     Szt(3,:)= 1i*kz(1,:,ic);                                                %Amplitude Xsi z +
    
    Uz(1,:) =-1i*kz(2,:,ic);                                                %Amplitude Phi sens z +
    Uz(2,:) = kr.^2;                                                        %Amplitude Psi z +
    Uz(3,:) = 0;                                                            %Amplitude Xsi z +
    
    Ur(1,:) = 1;                                                            %Amplitude Phi z +
    Ur(2,:) =-1i*kz(1,:,ic);                                                %Amplitude Psi z +
    Ur(3,:) = 1;                                                            %Amplitude Xsi z +
    
    %     Ut(1,:) =-1;                                                            %Amplitude Phi z +
    %     Ut(2,:) = 1i*kz(1,:,ic);                                                %Amplitude Psi z +
    %     Ut(3,:) =-1;                                                            %Amplitude Xsi z +
    
    %como gamma-=-gamma+ y nu-=-nu+, se deduce:
    % Szz(1,-)= Szz(1,+)
    % Szz(2,-)=-Szz(2,+)
    % Szz(3,-)= Szz(3,+)
    
    % Szr(1,-)=-Szr(1,+)
    % Szr(2,-)= Szr(2,+)
    % Szr(3,-)=-Szr(3,+)
    
    % Szt(1,-)=-Szt(1,+)
    % Szt(2,-)= Szt(2,+)
    % Szt(3,-)=-Szt(3,+)
    
    % Uz(1,-) =-Uz(1,+)
    % Uz(2,-) = Uz(2,+)
    % Uz(3,-) = Uz(3,+)
    
    % Ur(1,-) = Ur(1,+)
    % Ur(2,-) =-Ur(2,+)
    % Ur(3,-) = Ur(3,+)
    
    % Ut(1,-) = Ut(1,+)
    % Ut(2,-) =-Ut(2,+)
    % Ut(3,-) = Ut(3,+)
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMAzz EN 0
    A_DWN(1,1,:)= mic*Szz(1,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude Phi sens +
    A_DWN(1,2,:)= mic*Szz(2,:); %.*exp(-1i.*kz(2,:,ic).*0);             	%Amplitude Psi sens +
    
    %SIGMAzr EN 0
    A_DWN(2,1,:)= mic*Szr(1,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude Phi sens +
    A_DWN(2,2,:)= mic*Szr(2,:); %.*exp(-1i.*kz(2,:,ic).*0);                 %Amplitude Psi sens +
    
    B_DWN(1,1,:)= mic*Szr(3,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude Xsi sens +
    
    if ncapas>1
        %z=h des ondes descendantes avec origine en z=0
        exp1=exp(-1i.*kz(1,:,ic)*h);
        exp2=exp(-1i.*kz(2,:,ic)*h);
        %z=0 des ondes montantes avec origine en z=h
        exp3=exp1;%exp( 1i.*kz(1,:,ic).*(0-h));
        exp4=exp2;%exp( 1i.*kz(2,:,ic).*(0-h));
        
        %Amplitude des ondes avec kz sens -
        %SIGMAzz EN 0
        A_DWN(1,3,:)= mic*Szz(1,:).*exp4;                                   %Amplitude Phi sens -
        A_DWN(1,4,:)=-mic*Szz(2,:).*exp3;                                   %Amplitude Psi sens -
        
        %SIGMAzr EN 0
        A_DWN(2,3,:)=-mic*Szr(1,:).*exp4;                                   %Amplitude Phi sens -
        A_DWN(2,4,:)= mic*Szr(2,:).*exp3;                                   %Amplitude Psi sens -
        
        B_DWN(1,2,:)=-mic*Szr(3,:).*exp3;                                   %Amplitude Xsi sens -
        
        %Amplitude des ondes avec kz sens +
        %SIGMAzz EN h(ic)
        A_DWN(3,1,:)= mic*Szz(1,:).*exp2;                                   %Amplitude Phi sens +
        A_DWN(3,2,:)= mic*Szz(2,:).*exp1;                                   %Amplitude Psi sens +
        
        %SIGMAzr EN h(ic)
        A_DWN(4,1,:)= mic*Szr(1,:).*exp2;                                   %Amplitude Phi sens +
        A_DWN(4,2,:)= mic*Szr(2,:).*exp1;                                   %Amplitude Psi sens +
        
        B_DWN(2,1,:)= mic*Szr(3,:).*exp1;                                   %Amplitude Xsi sens +
        
        %Amplitude des ondes avec kz sens -
        %SIGMAzz EN h(ic)
        A_DWN(3,3,:)= mic*Szz(1,:);                                         %Amplitude Phi sens -
        A_DWN(3,4,:)=-mic*Szz(2,:);                                         %Amplitude Psi sens -
        
        %SIGMAzr EN h(ic)
        A_DWN(4,3,:)=-mic*Szr(1,:);                                         %Amplitude Phi sens -
        A_DWN(4,4,:)= mic*Szr(2,:);                                         %Amplitude Psi sens -
        
        B_DWN(2,2,:)=-mic*Szr(3,:);                                         %Amplitude Xsi sens -
        
        %Amplitude des ondes avec kz sens +
        %Uz EN h(ic)
        A_DWN(5,1,:)= mic*Uz(1,:).*exp2;                                    %Amplitude Phi sens +
        A_DWN(5,2,:)= mic*Uz(2,:).*exp1;                                    %Amplitude Psi sens +
        
        %Ur EN h(ic)
        A_DWN(6,1,:)= mic*Ur(1,:).*exp2;                                    %Amplitude Phi sens +
        A_DWN(6,2,:)= mic*Ur(2,:).*exp1;                                    %Amplitude Psi sens +
        
        B_DWN(3,1,:)= mic*Ur(3,:).*exp1;                                    %Amplitude Xsi sens +
        
        %Amplitude des ondes avec kz sens -
        %Uz EN h(ic)
        A_DWN(5,3,:)=-mic*Uz(1,:);                                          %Amplitude Phi sens -
        A_DWN(5,4,:)= mic*Uz(2,:);                                          %Amplitude Psi sens -
        
        %Ur EN h(ic)
        A_DWN(6,3,:)= mic*Ur(1,:);                                          %Amplitude Phi sens -
        A_DWN(6,4,:)=-mic*Ur(2,:);                                          %Amplitude Psi sens -
        
        B_DWN(3,2,:)= mic*Ur(3,:);                                          %Amplitude Xsi sens -
    end
end

for ic=2:(ncapas-1)
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    h  = para.reg(1).sub(ic).h;
    xi = kz(1,:,ic).^2-kr.^2;
    
    Szz(1,:)=-xi;                                                           %Amplitude Phi sens z +
    Szz(2,:)=-2i*kz(1,:,ic).*kr.^2;                                          %Amplitude Psi z +
    Szz(3,:)= 0;                                                            %Amplitude Xsi z +
    
    Szr(1,:)= 2i*kz(2,:,ic);                                                %Amplitude Phi z +
    Szr(2,:)= xi;                                                           %Amplitude Psi z +
    Szr(3,:)= 1i*kz(1,:,ic);                                                %Amplitude Xsi z +
    
    %     Szt(1,:)= 2i*kz(2,:,ic);                                                %Amplitude Phi z +
    %     Szt(2,:)= xi;                                                           %Amplitude Psi z +
    %     Szt(3,:)= 1i*kz(1,:,ic);                                                %Amplitude Xsi z +
    
    Szz     = MAT(ic).Ci(6,6)*Szz;
    Szr     = MAT(ic).Ci(6,6)*Szr;
    
    Uz(1,:) =-1i*kz(2,:,ic);                                                %Amplitude Phi sens z +
    Uz(2,:) = kr.^2;                                                        %Amplitude Psi z +
    Uz(3,:) = 0;                                                            %Amplitude Xsi z +
    
    Ur(1,:) = 1;                                                            %Amplitude Phi z +
    Ur(2,:) =-1i*kz(1,:,ic);                                                %Amplitude Psi z +
    Ur(3,:) = 1;                                                            %Amplitude Xsi z +
    
    %     Ut(1,:) = 1;                                                            %Amplitude Phi z +
    %     Ut(2,:) =-1i*kz(1,:,ic);                                                %Amplitude Psi z +
    %     Ut(3,:) = 1;                                                            %Amplitude Xsi z +
    
    %z=h des ondes descendantes avec origine en z=0
    exp1=exp(-1i*kz(1,:,ic)*h);
    exp2=exp(-1i*kz(2,:,ic)*h);
    %z=0 des ondes montantes avec origine en z=h
    exp3=exp1;%exp( 1i.*kz(1,:,ic).*(0-h));
    exp4=exp2;%exp( 1i.*kz(2,:,ic).*(0-h));
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMAzz EN 0
    A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)= mic*Szz(1,:);                         %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+1,(ic-1)*4+2,:)= mic*Szz(2,:);                         %Amplitude Psi sens +
    
    %SIGMAzr EN 0
    A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)= mic*Szr(1,:);                         %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+2,(ic-1)*4+2,:)= mic*Szr(2,:);                         %Amplitude Psi sens +
    
    B_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*Szr(3,:);                         %Amplitude Xsi sens +
    
    %Amplitude des ondes avec kz sens -
    %SIGMAzz EN 0
    A_DWN((ic-2)*4+2+1,(ic-1)*4+3,:)= mic*Szz(1,:).*exp4;                   %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+1,(ic-1)*4+4,:)=-mic*Szz(2,:).*exp3;                   %Amplitude Psi sens -
    
    %SIGMAzr EN 0
    A_DWN((ic-2)*4+2+2,(ic-1)*4+3,:)=-mic*Szr(1,:).*exp4;                   %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+2,(ic-1)*4+4,:)= mic*Szr(2,:).*exp3;                   %Amplitude Psi sens -
    
    B_DWN((ic-2)*2+1+1,(ic-1)*2+2,:)=-mic*Szr(3,:).*exp3;                   %Amplitude Xsi sens -
    
    %Amplitude des ondes avec kz sens +
    %Uz EN 0(ic)
    A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)= mic*Uz(1,:);                          %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+3,(ic-1)*4+2,:)= mic*Uz(2,:);                          %Amplitude Psi sens +
    
    %Ur EN 0(ic)
    A_DWN((ic-2)*4+2+4,(ic-1)*4+1,:)= mic*Ur(1,:);                          %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+4,(ic-1)*4+2,:)= mic*Ur(2,:);                          %Amplitude Psi sens +
    
    B_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic*Ur(3,:);                          %Amplitude Xsi sens +
    
    %Amplitude des ondes avec kz sens -
    %Uz EN 0(ic)
    A_DWN((ic-2)*4+2+3,(ic-1)*4+3,:)=-mic*Uz(1,:).*exp4;                    %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+3,(ic-1)*4+4,:)= mic*Uz(2,:).*exp3;                    %Amplitude Psi sens -
    
    %Ur EN 0(ic)
    A_DWN((ic-2)*4+2+4,(ic-1)*4+3,:)= mic*Ur(1,:).*exp4;                    %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+4,(ic-1)*4+4,:)=-mic*Ur(2,:).*exp3;                    %Amplitude Psi sens -
    
    B_DWN((ic-2)*2+1+2,(ic-1)*2+2,:)= mic*Ur(3,:).*exp3;                    %Amplitude Xsi sens -
    
    %Amplitude des ondes avec kz sens +
    %SIGMAzz EN h(ic)
    A_DWN((ic-2)*4+2+5,(ic-1)*4+1,:)= mic*Szz(1,:).*exp2;                   %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+5,(ic-1)*4+2,:)= mic*Szz(2,:).*exp1;                   %Amplitude Psi sens +
    
    %SIGMAzr EN h(ic)
    A_DWN((ic-2)*4+2+6,(ic-1)*4+1,:)= mic*Szr(1,:).*exp2;                   %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+6,(ic-1)*4+2,:)= mic*Szr(2,:).*exp1;                   %Amplitude Psi sens +
    
    B_DWN((ic-2)*2+1+3,(ic-1)*2+1,:)= mic*Szr(3,:).*exp1;                   %Amplitude Xsi sens +
    
    %Amplitude des ondes avec kz sens -
    %SIGMAzz EN h(ic)
    A_DWN((ic-2)*4+2+5,(ic-1)*4+3,:)= mic*Szz(1,:);                         %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+5,(ic-1)*4+4,:)=-mic*Szz(2,:);                         %Amplitude Psi sens -
    
    %SIGMAzr EN h(ic)
    A_DWN((ic-2)*4+2+6,(ic-1)*4+3,:)=-mic*Szr(1,:);                         %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+6,(ic-1)*4+4,:)= mic*Szr(2,:);                         %Amplitude Psi sens -
    
    B_DWN((ic-2)*2+1+3,(ic-1)*2+2,:)=-mic*Szr(3,:);                         %Amplitude Xsi sens -
    
    %Amplitude des ondes avec kz sens +
    %Uz EN h(ic)
    A_DWN((ic-2)*4+2+7,(ic-1)*4+1,:)= mic*Uz(1,:).*exp2;                   %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+7,(ic-1)*4+2,:)= mic*Uz(2,:).*exp1;                   %Amplitude Psi sens +
    
    %Ur EN h(ic)
    A_DWN((ic-2)*4+2+8,(ic-1)*4+1,:)= mic*Ur(1,:).*exp2;                   %Amplitude Phi sens +
    A_DWN((ic-2)*4+2+8,(ic-1)*4+2,:)= mic*Ur(2,:).*exp1;                   %Amplitude Psi sens +
    
    B_DWN((ic-2)*2+1+4,(ic-1)*2+1,:)= mic*Ur(3,:).*exp1;                    %Amplitude Xsi sens +
    
    %Amplitude des ondes avec kz sens -
    %Uz EN h(ic)
    A_DWN((ic-2)*4+2+7,(ic-1)*4+3,:)=-mic*Uz(1,:);                         %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+7,(ic-1)*4+4,:)= mic*Uz(2,:);                         %Amplitude Psi sens -
    
    %Ur EN h(ic)
    A_DWN((ic-2)*4+2+8,(ic-1)*4+3,:)= mic*Ur(1,:);                         %Amplitude Phi sens -
    A_DWN((ic-2)*4+2+8,(ic-1)*4+4,:)=-mic*Ur(2,:);                         %Amplitude Psi sens -
    
    B_DWN((ic-2)*2+1+4,(ic-1)*2+2,:)= mic*Ur(3,:);                          %Amplitude Xsi sens -
end

if ncapas>1
    for ic=ncapas
        %***********************************************************%
        %   calcul de la contrainte issue de la solution homogene   %
        %***********************************************************%
        xi = kz(1,:,ic).^2-kr.^2;
        
        Szz(1,:)=-xi;                                                       %Amplitude Phi sens z +
        Szz(2,:)=-2i*kz(1,:,ic).*kr.^2;                                  	%Amplitude Psi z +
        Szz(3,:)= 0;                                                    	%Amplitude Xsi z +
        
        Szr(1,:)= 2i*kz(2,:,ic);                                         	%Amplitude Phi z +
        Szr(2,:)= xi;                                                   	%Amplitude Psi z +
        Szr(3,:)= 1i*kz(1,:,ic);                                        	%Amplitude Xsi z +
        
        %         Szt(1,:)= 2i*kz(2,:,ic);                                            %Amplitude Phi z +
        %         Szt(2,:)= xi;                                                   	%Amplitude Psi z +
        %         Szt(3,:)= 1i*kz(1,:,ic);                                        	%Amplitude Xsi z +
        
        Szz     = MAT(ic).Ci(6,6)*Szz;
        Szr     = MAT(ic).Ci(6,6)*Szr;
        
        
        Uz(1,:) =-1i*kz(2,:,ic);                                         	%Amplitude Phi sens z +
        Uz(2,:) = kr.^2;                                                  	%Amplitude Psi z +
        Uz(3,:) = 0;                                                     	%Amplitude Xsi z +
        
        Ur(1,:) = 1;                                                      	%Amplitude Phi z +
        Ur(2,:) =-1i*kz(1,:,ic);                                        	%Amplitude Psi z +
        Ur(3,:) = 1;                                                     	%Amplitude Xsi z +
        
        %         Ut(1,:) = 1;                                                     	%Amplitude Phi z +
        %         Ut(2,:) =-1i*kz(1,:,ic);                                          	%Amplitude Psi z +
        %         Ut(3,:) = 1;                                                      	%Amplitude Xsi z +
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        mic	= (-1)^ic;
        
        %SIGMAzz EN 0
        A_DWN((ic-2)*4+2+1,(ic-1)*4+1,:)= mic*Szz(1,:);                     %Amplitude Phi sens +
        A_DWN((ic-2)*4+2+1,(ic-1)*4+2,:)= mic*Szz(2,:);                     %Amplitude Psi sens +
        
        %SIGMAzr EN 0
        A_DWN((ic-2)*4+2+2,(ic-1)*4+1,:)= mic*Szr(1,:);                     %Amplitude Phi sens +
        A_DWN((ic-2)*4+2+2,(ic-1)*4+2,:)= mic*Szr(2,:);                     %Amplitude Psi sens +
        
        B_DWN((ic-2)*2+1+1,(ic-1)*2+1,:)= mic*Szr(3,:);                     %Amplitude Xsi sens +
        
        %Amplitude des ondes avec kz sens +
        %Uz EN 0(ic)
        A_DWN((ic-2)*4+2+3,(ic-1)*4+1,:)= mic*Uz(1,:);                      %Amplitude Phi sens +
        A_DWN((ic-2)*4+2+3,(ic-1)*4+2,:)= mic*Uz(2,:);                      %Amplitude Psi sens +
        
        %Ur EN 0(ic)
        A_DWN((ic-2)*4+2+4,(ic-1)*4+1,:)= mic*Ur(1,:);                      %Amplitude Phi sens +
        A_DWN((ic-2)*4+2+4,(ic-1)*4+2,:)= mic*Ur(2,:);                      %Amplitude Psi sens +
        
        B_DWN((ic-2)*2+1+2,(ic-1)*2+1,:)= mic*Ur(3,:);                      %Amplitude Xsi sens +
    end
end
% end % spmd

DWN.A_DWN = A_DWN;%gather(A_DWN);
DWN.B_DWN = B_DWN;%gather(B_DWN);
% toc
% disp(A_DWN(:,:,round(end/2)))
% error('debugging')
end