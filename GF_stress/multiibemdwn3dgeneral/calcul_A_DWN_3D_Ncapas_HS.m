function DWN=calcul_A_DWN_3D_Ncapas_HS(para,DWN)
% Matrices A,B del medio estratificado (DWN) solución en cartesianas
% (lento)
% METHODE : @ da los valores de kz (componente vertical de P,SV,SH) por un medio isotropo
%           @ da los vectores propio, polarizacion associados a kx, ky y kz
%           @ calculo de la matrice de los coefficiente de la sol homogenea
%             a partir de las condiciones en las interfases
% OJO!!! se cambia los ejes con respectos a los modelos PSV 1=x, 2==y, 3==z

ncapas      = para.nsubmed;
kx          = DWN.kx;
ky          = DWN.kx;


kr          = sqrt(kx.^2+ky.^2);
nk2       	= length(kx(1,:));
if nk2==1
    nk2=length(DWN.omegac);
end
u1          = zeros(3,nk2,nk2,ncapas);
u2          = zeros(3,nk2,nk2,ncapas);
u3          = zeros(3,nk2,nk2,ncapas);
kz          = zeros(3,nk2,nk2,ncapas);

S31         = zeros(3,nk2,nk2);
S32         = zeros(3,nk2,nk2);
S33         = zeros(3,nk2,nk2);

A_DWN       = zeros(6*(para.nsubmed-1)+3,6*(para.nsubmed-1)+3,nk2,nk2);

%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%
for ic=1:ncapas
    ksi         = para.reg(1).sub(ic).ksi;
    kpi         = para.reg(1).sub(ic).kpi;
    kz(1,:,ic)  = sqrt(ksi.^2-kx.^2-ky.'.^2); %k1S
    kz(2,:,ic) 	= sqrt(kpi.^2-kx.^2-ky.'.^2); %k1P
    kz(3,:,ic)  = kz(1,:,ic);
end
kz  =-kz.*(sign(imag(kz)).*(imag(kz)~=0)+(imag(kz)==0));

for ic=1:ncapas
    %polarizacion normalizadas
    u1(1,:,ic)  = kz(1,:,ic).*kx./kr/ksi;	% SV ux
    u2(1,:,ic)  = kz(1,:,ic).*ky./kr/ksi;	% SV uy
    u3(1,:,ic)  =-kr/ksi;                   % SV uz
    
    u1(2,:,ic)  = kx /kpi;                  % P ux
    u2(2,:,ic)  = ky /kpi;                  % P uy
    u3(2,:,ic)  = kz(2,:,ic)/kpi;       	% P uz
    
    u1(3,:,ic)  = ky./kr;                   % SH ux
    u2(3,:,ic)  =-kx./kr;                   % SH uy
    u3(3,:,ic)  = 0;                        % SH uz
end
DWN.kz=kz;
DWN.u1=u1;
DWN.u2=u2;
DWN.u3=u3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%les ondes en exp(-1i.*kz(:,ic).*z    ) ont leur origine en z=0 interface supérieure
%les ondes en exp( 1i.*kz(:,ic).*(z-h)) ont leur origine en z=h interface inférieure

for ic=1
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    Ci = para.reg(1).sub(ic).Ci;
    h  = para.reg(1).sub(ic).h;
    for i=1:3
        % S31(i,:)=-1i*Ci(5,5)*(kz(i,:,ic).*u1(i,:,ic) + kx.*u3(i,:,ic));                                  	%Amplitude longi & transverse sens z +
        % S32(i,:)=-1i*Ci(4,4)*(kz(i,:,ic).*u2(i,:,ic) + ky.*u3(i,:,ic));                                  	%Amplitude longi & transverse sens z +
        % S33(i,:)=-1i*(kx(i,:,ic).*Ci(3,1).*u1(i,:,ic) + ky.*Ci(3,2).*u2(i,:,ic)+kz.*Ci(3,3).*u3(i,:,ic));	%Amplitude longi & transverse sens z +
        % reescritura tomando en cuenta las simetrias de Cij
        S31(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u1(i,:,ic) + kx.*u3(i,:,ic));
        S32(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u2(i,:,ic) + ky.*u3(i,:,ic));
        S33(i,:)=-1i*(Ci(1,2)*(kx.*u1(i,:,ic) + ky.*u2(i,:,ic))+Ci(1,1)*kz(i,:,ic).*u3(i,:,ic));
    end
    %se deduce que:
    % S31(1,-)= S31(1,+) pour kzSV-= kzSV+
    % S32(1,-)= S32(1,+)
    % S33(1,-)=-S33(1,+)
    
    % S31(2,-)=-S31(2,+) pour kz(2,:,ic)- = kz(2,:,ic)+
    % S32(2,-)=-S32(2,+)
    % S33(2,-)= S33(2,+)
    
    % S31(3,-)=-S31(3,+) pour kzSH-= kzSH+
    % S32(3,-)=-S32(3,+)
    % S33(3,-)= S33(3,+)
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMA31 EN 0
    A_DWN(1,1,:)=mic*S31(1,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    A_DWN(1,2,:)=mic*S31(2,:); %.*exp(-1i.*kz(2,:,ic).*0);                 %Amplitude longi sens +
    A_DWN(1,3,:)=mic*S31(3,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    
    %SIGMA32 EN 0
    A_DWN(2,1,:)=mic*S32(1,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    A_DWN(2,2,:)=mic*S32(2,:); %.*exp(-1i.*kz(2,:,ic).*0);                 %Amplitude longi sens +
    A_DWN(2,3,:)=mic*S32(3,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    
    %SIGMA33 EN 0
    A_DWN(3,1,:)=mic*S33(1,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    A_DWN(3,2,:)=mic*S33(2,:); %.*exp(-1i.*kz(2,:,ic).*0);                 %Amplitude longi sens +
    A_DWN(3,3,:)=mic*S33(3,:); %.*exp(-1i.*kz(1,:,ic).*0);                 %Amplitude trans sens +
    
    if ncapas>1
        %z=h des ondes descendantes avec origine en z=0
        exp1=exp(-1i.*kz(1,:,ic).*h);
        exp2=exp(-1i.*kz(2,:,ic).*h);
        %z=0 des ondes montantes avec origine en z=h
        exp3=exp1;%exp( 1i.*kz(1,:,ic).*(0-h));
        exp4=exp2;%exp( 1i.*kz(2,:,ic).*(0-h));
        
        %Amplitude des ondes avec kz sens -
        %SIGMA31 EN 0
        A_DWN(1,4,:)= mic*S31(1,:).*exp3;                                  %Amplitude trans sens -
        A_DWN(1,5,:)=-mic*S31(2,:).*exp4;                                  %Amplitude longi sens -
        A_DWN(1,6,:)=-mic*S31(3,:).*exp3;                                  %Amplitude trans sens -
        
        %SIGMA32 EN 0
        A_DWN(2,4,:)= mic*S32(1,:).*exp3;                                  %Amplitude trans sens -
        A_DWN(2,5,:)=-mic*S32(2,:).*exp4;                                  %Amplitude longi sens -
        A_DWN(2,6,:)=-mic*S32(3,:).*exp3;                                  %Amplitude trans sens -
        
        %SIGMA33 EN 0
        A_DWN(3,4,:)=-mic*S33(1,:).*exp3;                                  %Amplitude trans sens -
        A_DWN(3,5,:)= mic*S33(2,:).*exp4;                                  %Amplitude longi sens -
        A_DWN(3,6,:)= mic*S33(3,:).*exp3;                                  %Amplitude trans sens -
        
        %Amplitude des ondes avec kz sens +
        %SIGMA31 EN h(ic)
        A_DWN(4,1,:)= mic*S31(1,:).*exp1;                                  %Amplitude trans sens +
        A_DWN(4,2,:)= mic*S31(2,:).*exp2;                                  %Amplitude longi sens +
        A_DWN(4,3,:)= mic*S31(3,:).*exp1;                                  %Amplitude trans sens +
        
        %SIGMA32 EN h(ic)
        A_DWN(5,1,:)= mic*S32(1,:).*exp1;                                  %Amplitude trans sens +
        A_DWN(5,2,:)= mic*S32(2,:).*exp2;                                  %Amplitude longi sens +
        A_DWN(5,3,:)= mic*S32(3,:).*exp1;                                  %Amplitude trans sens +
        
        %SIGMA33 EN h(ic)
        A_DWN(6,1,:)= mic*S33(1,:).*exp1;                                  %Amplitude trans sens +
        A_DWN(6,2,:)= mic*S33(2,:).*exp2;                                  %Amplitude longi sens +
        A_DWN(6,3,:)= mic*S33(3,:).*exp1;                                  %Amplitude trans sens +
        
        %Amplitude des ondes avec kz sens -
        %SIGMA31 EN h(ic)
        A_DWN(4,4,:)= mic*S31(1,:);                                        %Amplitude trans sens -
        A_DWN(4,5,:)=-mic*S31(2,:);                                        %Amplitude longi sens -
        A_DWN(4,6,:)=-mic*S31(3,:);                                        %Amplitude trans sens -
        
        %SIGMA32 EN h(ic)
        A_DWN(5,4,:)= mic*S32(1,:);                                        %Amplitude trans sens -
        A_DWN(5,5,:)=-mic*S32(2,:);                                        %Amplitude longi sens -
        A_DWN(5,6,:)=-mic*S32(3,:);                                        %Amplitude trans sens -
        
        %SIGMA33 EN h(ic)
        A_DWN(6,4,:)=-mic*S33(1,:);                                        %Amplitude trans sens -
        A_DWN(6,5,:)= mic*S33(2,:);                                        %Amplitude longi sens -
        A_DWN(6,6,:)= mic*S33(3,:);                                        %Amplitude trans sens -
        
        %Amplitude des ondes avec kz sens +
        %U1 EN h(ic)
        A_DWN(7,1,:)= mic*u1(1,:,ic).*exp1;                                %Amplitude trans sens +
        A_DWN(7,2,:)= mic*u1(2,:,ic).*exp2;                                %Amplitude longi sens +
        A_DWN(7,3,:)= mic*u1(3,:,ic).*exp1;                                %Amplitude trans sens +
        
        %U2 EN h(ic)
        A_DWN(8,1,:)= mic*u2(1,:,ic).*exp1;                                %Amplitude trans sens +
        A_DWN(8,2,:)= mic*u2(2,:,ic).*exp2;                                %Amplitude longi sens +
        A_DWN(8,3,:)= mic*u2(3,:,ic).*exp1;                                %Amplitude trans sens +
        
        %U3 EN h(ic)
        A_DWN(9,1,:)= mic*u3(1,:,ic).*exp1;                                %Amplitude trans sens +
        A_DWN(9,2,:)= mic*u3(2,:,ic).*exp2;                                %Amplitude longi sens +
        A_DWN(9,3,:)= mic*u3(3,:,ic).*exp1;                                %Amplitude trans sens +
        
        %Amplitude des ondes avec kz sens -
        %U1 EN h(ic)
        A_DWN(7,4,:)=-mic*u1(1,:,ic);                                      %Amplitude trans sens -
        A_DWN(7,5,:)= mic*u1(2,:,ic);                                      %Amplitude longi sens -
        A_DWN(7,6,:)= mic*u1(3,:,ic);                                      %Amplitude trans sens -
        
        %U2 EN h(ic)
        A_DWN(8,4,:)=-mic*u2(1,:,ic);                                      %Amplitude trans sens -
        A_DWN(8,5,:)= mic*u2(2,:,ic);                                      %Amplitude longi sens -
        A_DWN(8,6,:)= mic*u2(3,:,ic);                                      %Amplitude trans sens -
        
        %U3 EN h(ic)
        A_DWN(9,4,:)= mic*u3(1,:,ic);                                      %Amplitude trans sens -
        A_DWN(9,5,:)=-mic*u3(2,:,ic);                                      %Amplitude longi sens -
        A_DWN(9,6,:)= mic*u3(3,:,ic);                                      %Amplitude trans sens -
    end
end

for ic=2:(ncapas-1)
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    Ci = para.reg(1).sub(ic).Ci;
    h  = para.reg(1).sub(ic).h;
    
    for i=1:3
        S31(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u1(i,:,ic) + kx.*u3(i,:,ic));
        S32(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u2(i,:,ic) + ky.*u3(i,:,ic));
        S33(i,:)=-1i*(Ci(1,2)*(kx.*u1(i,:,ic) + ky.*u2(i,:,ic))+Ci(1,1)*kz(i,:,ic).*u3(i,:,ic));
    end
    
    %z=h des ondes descendantes avec origine en z=0
    exp1=exp(-1i.*kz(1,:,ic).*h);
    exp2=exp(-1i.*kz(2,:,ic).*h);
    %z=0 des ondes montantes avec origine en z=h
    exp3=exp1;%exp( 1i.*kz(1,:,ic).*(0-h));
    exp4=exp2;%exp( 1i.*kz(2,:,ic).*(0-h));
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMA31 EN 0
    A_DWN((ic-2)*6+3+1,(ic-1)*6+1,:)=mic*S31(1,:);                         %Amplitude trans sens +
    A_DWN((ic-2)*6+3+1,(ic-1)*6+2,:)=mic*S31(2,:);                         %Amplitude longi sens +
    A_DWN((ic-2)*6+3+1,(ic-1)*6+3,:)=mic*S31(3,:);                         %Amplitude trans sens +
    
    %SIGMA32 EN 0
    A_DWN((ic-2)*6+3+2,(ic-1)*6+1,:)=mic*S32(1,:);                         %Amplitude trans sens +
    A_DWN((ic-2)*6+3+2,(ic-1)*6+2,:)=mic*S32(2,:);                         %Amplitude longi sens +
    A_DWN((ic-2)*6+3+2,(ic-1)*6+3,:)=mic*S32(3,:);                         %Amplitude trans sens +
    
    %SIGMA33 EN 0
    A_DWN((ic-2)*6+3+3,(ic-1)*6+1,:)=mic*S33(1,:);                         %Amplitude trans sens +
    A_DWN((ic-2)*6+3+3,(ic-1)*6+2,:)=mic*S33(2,:);                         %Amplitude longi sens +
    A_DWN((ic-2)*6+3+3,(ic-1)*6+3,:)=mic*S33(3,:);                         %Amplitude trans sens +
    
    
    %Amplitude des ondes avec kz sens -
    %SIGMA31 EN 0
    A_DWN((ic-2)*6+3+1,(ic-1)*6+4,:)= mic*S31(1,:).*exp3;                  %Amplitude trans sens -
    A_DWN((ic-2)*6+3+1,(ic-1)*6+5,:)=-mic*S31(2,:).*exp4;                  %Amplitude longi sens -
    A_DWN((ic-2)*6+3+1,(ic-1)*6+6,:)=-mic*S31(3,:).*exp3;                  %Amplitude trans sens -
    
    %SIGMA32 EN 0
    A_DWN((ic-2)*6+3+2,(ic-1)*6+4,:)= mic*S32(1,:).*exp3;                  %Amplitude trans sens -
    A_DWN((ic-2)*6+3+2,(ic-1)*6+5,:)=-mic*S32(2,:).*exp4;                  %Amplitude longi sens -
    A_DWN((ic-2)*6+3+2,(ic-1)*6+6,:)=-mic*S32(3,:).*exp3;                  %Amplitude trans sens -
    
    %SIGMA33 EN 0
    A_DWN((ic-2)*6+3+3,(ic-1)*6+4,:)=-mic*S33(1,:).*exp3;                  %Amplitude trans sens -
    A_DWN((ic-2)*6+3+3,(ic-1)*6+5,:)= mic*S33(2,:).*exp4;                  %Amplitude longi sens -
    A_DWN((ic-2)*6+3+3,(ic-1)*6+6,:)= mic*S33(3,:).*exp3;                  %Amplitude trans sens -
    
    %Amplitude des ondes avec kz sens +
    %U1 EN 0(ic)
    A_DWN((ic-2)*6+3+4,(ic-1)*6+1,:)= mic*u1(1,:,ic);                      %Amplitude trans sens +
    A_DWN((ic-2)*6+3+4,(ic-1)*6+2,:)= mic*u1(2,:,ic);                      %Amplitude longi sens +
    A_DWN((ic-2)*6+3+4,(ic-1)*6+3,:)= mic*u1(3,:,ic);                      %Amplitude trans sens +
    
    %U2 EN 0(ic)
    A_DWN((ic-2)*6+3+5,(ic-1)*6+1,:)= mic*u2(1,:,ic);                      %Amplitude trans sens +
    A_DWN((ic-2)*6+3+5,(ic-1)*6+2,:)= mic*u2(2,:,ic);                      %Amplitude longi sens +
    A_DWN((ic-2)*6+3+5,(ic-1)*6+3,:)= mic*u2(3,:,ic);                      %Amplitude trans sens +
    
    %U3 EN 0(ic)
    A_DWN((ic-2)*6+3+6,(ic-1)*6+1,:)= mic*u3(1,:,ic);                      %Amplitude trans sens +
    A_DWN((ic-2)*6+3+6,(ic-1)*6+2,:)= mic*u3(2,:,ic);                      %Amplitude longi sens +
    A_DWN((ic-2)*6+3+6,(ic-1)*6+3,:)= mic*u3(3,:,ic);                      %Amplitude trans sens +
    
    %Amplitude des ondes avec kz sens -
    %U1 EN 0(ic)
    A_DWN((ic-2)*6+3+4,(ic-1)*6+4,:)=-mic*u1(1,:,ic).*exp3;                %Amplitude trans sens -
    A_DWN((ic-2)*6+3+4,(ic-1)*6+5,:)= mic*u1(2,:,ic).*exp4;                %Amplitude longi sens -
    A_DWN((ic-2)*6+3+4,(ic-1)*6+6,:)= mic*u1(3,:,ic).*exp3;                %Amplitude trans sens -
    
    %U2 EN 0(ic)
    A_DWN((ic-2)*6+3+5,(ic-1)*6+4,:)=-mic*u2(1,:,ic).*exp3;                %Amplitude trans sens -
    A_DWN((ic-2)*6+3+5,(ic-1)*6+5,:)= mic*u2(2,:,ic).*exp4;                %Amplitude longi sens -
    A_DWN((ic-2)*6+3+5,(ic-1)*6+6,:)= mic*u2(3,:,ic).*exp3;                %Amplitude trans sens -
    
    %U3 EN 0(ic)
    A_DWN((ic-2)*6+3+6,(ic-1)*6+4,:)= mic*u3(1,:,ic).*exp3;                %Amplitude trans sens -
    A_DWN((ic-2)*6+3+6,(ic-1)*6+5,:)=-mic*u3(2,:,ic).*exp4;                %Amplitude longi sens -
    A_DWN((ic-2)*6+3+6,(ic-1)*6+6,:)= mic*u3(3,:,ic).*exp3;                %Amplitude trans sens -
    
    %Amplitude des ondes avec kz sens +
    %SIGMA31 EN h(ic)
    A_DWN((ic-2)*6+3+7,(ic-1)*6+1,:)= mic*S31(1,:).*exp1;                  %Amplitude trans sens +
    A_DWN((ic-2)*6+3+7,(ic-1)*6+2,:)= mic*S31(2,:).*exp2;                  %Amplitude longi sens +
    A_DWN((ic-2)*6+3+7,(ic-1)*6+3,:)= mic*S31(3,:).*exp1;                  %Amplitude trans sens +
    
    %SIGMA32 EN h(ic)
    A_DWN((ic-2)*6+3+8,(ic-1)*6+1,:)= mic*S32(1,:).*exp1;                  %Amplitude trans sens +
    A_DWN((ic-2)*6+3+8,(ic-1)*6+2,:)= mic*S32(2,:).*exp2;                  %Amplitude longi sens +
    A_DWN((ic-2)*6+3+8,(ic-1)*6+3,:)= mic*S32(3,:).*exp1;                  %Amplitude trans sens +
    
    %SIGMA33 EN h(ic)
    A_DWN((ic-2)*6+3+9,(ic-1)*6+1,:)= mic*S33(1,:).*exp1;                  %Amplitude trans sens +
    A_DWN((ic-2)*6+3+9,(ic-1)*6+2,:)= mic*S33(2,:).*exp2;                  %Amplitude longi sens +
    A_DWN((ic-2)*6+3+9,(ic-1)*6+3,:)= mic*S33(3,:).*exp1;                  %Amplitude trans sens +
    
    %Amplitude des ondes avec kz sens -
    %SIGMA31 EN h(ic)
    A_DWN((ic-2)*6+3+7,(ic-1)*6+4,:)= mic*S31(1,:);                        %Amplitude trans sens -
    A_DWN((ic-2)*6+3+7,(ic-1)*6+5,:)=-mic*S31(2,:);                        %Amplitude longi sens -
    A_DWN((ic-2)*6+3+7,(ic-1)*6+6,:)=-mic*S31(3,:);                        %Amplitude trans sens -
    
    %SIGMA32 EN h(ic)
    A_DWN((ic-2)*6+3+8,(ic-1)*6+4,:)= mic*S32(1,:);                        %Amplitude trans sens -
    A_DWN((ic-2)*6+3+8,(ic-1)*6+5,:)=-mic*S32(2,:);                        %Amplitude longi sens -
    A_DWN((ic-2)*6+3+8,(ic-1)*6+6,:)=-mic*S32(3,:);                        %Amplitude trans sens -
    
    %SIGMA33 EN h(ic)
    A_DWN((ic-2)*6+3+9,(ic-1)*6+4,:)=-mic*S33(1,:);                        %Amplitude trans sens -
    A_DWN((ic-2)*6+3+9,(ic-1)*6+5,:)= mic*S33(2,:);                        %Amplitude longi sens -
    A_DWN((ic-2)*6+3+9,(ic-1)*6+6,:)= mic*S33(3,:);                        %Amplitude trans sens -
    
    %Amplitude des ondes avec kz sens +
    %U1 EN h(ic)
    A_DWN((ic-2)*6+3+10,(ic-1)*6+1,:)= mic*u1(1,:,ic).*exp1;               %Amplitude trans sens +
    A_DWN((ic-2)*6+3+10,(ic-1)*6+2,:)= mic*u1(2,:,ic).*exp2;               %Amplitude longi sens +
    A_DWN((ic-2)*6+3+10,(ic-1)*6+3,:)= mic*u1(3,:,ic).*exp1;               %Amplitude trans sens +
    
    %U2 EN h(ic)
    A_DWN((ic-2)*6+3+11,(ic-1)*6+1,:)= mic*u2(1,:,ic).*exp1;               %Amplitude trans sens +
    A_DWN((ic-2)*6+3+11,(ic-1)*6+2,:)= mic*u2(2,:,ic).*exp2;               %Amplitude longi sens +
    A_DWN((ic-2)*6+3+11,(ic-1)*6+3,:)= mic*u2(3,:,ic).*exp1;               %Amplitude trans sens +
    
    %U3 EN h(ic)
    A_DWN((ic-2)*6+3+12,(ic-1)*6+1,:)= mic*u3(1,:,ic).*exp1;               %Amplitude trans sens +
    A_DWN((ic-2)*6+3+12,(ic-1)*6+2,:)= mic*u3(2,:,ic).*exp2;               %Amplitude longi sens +
    A_DWN((ic-2)*6+3+12,(ic-1)*6+3,:)= mic*u3(3,:,ic).*exp1;               %Amplitude trans sens +
    
    %Amplitude des ondes avec kz sens -
    %U1 EN h(ic)
    A_DWN((ic-2)*6+3+10,(ic-1)*6+4,:)=-mic*u1(1,:,ic);                     %Amplitude trans sens -
    A_DWN((ic-2)*6+3+10,(ic-1)*6+5,:)= mic*u1(2,:,ic);                     %Amplitude longi sens -
    A_DWN((ic-2)*6+3+10,(ic-1)*6+6,:)= mic*u1(3,:,ic);                     %Amplitude trans sens -
    
    %U2 EN h(ic)
    A_DWN((ic-2)*6+3+11,(ic-1)*6+4,:)=-mic*u2(1,:,ic);                     %Amplitude trans sens -
    A_DWN((ic-2)*6+3+11,(ic-1)*6+5,:)= mic*u2(2,:,ic);                     %Amplitude longi sens -
    A_DWN((ic-2)*6+3+11,(ic-1)*6+6,:)= mic*u2(3,:,ic);                     %Amplitude trans sens -
    
    %U3 EN h(ic)
    A_DWN((ic-2)*6+3+12,(ic-1)*6+4,:)= mic*u3(1,:,ic);                     %Amplitude trans sens -
    A_DWN((ic-2)*6+3+12,(ic-1)*6+5,:)=-mic*u3(2,:,ic);                     %Amplitude longi sens -
    A_DWN((ic-2)*6+3+12,(ic-1)*6+6,:)= mic*u3(3,:,ic);                     %Amplitude trans sens -
end

if ncapas>1
    for ic=ncapas
        %***********************************************************%
        %   calcul de la contrainte issue de la solution homogene   %
        %***********************************************************%
        Ci = para.reg(1).sub(ic).Ci;
        
        for i=1:3
            S31(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u1(i,:,ic) + kx.*u3(i,:,ic));
            S32(i,:)=-1i*Ci(6,6)*(kz(i,:,ic).*u2(i,:,ic) + ky.*u3(i,:,ic));
            S33(i,:)=-1i*(Ci(1,2)*(kx.*u1(i,:,ic) + ky.*u2(i,:,ic))+Ci(1,1)*kz(i,:,ic).*u3(i,:,ic));
        end
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        mic	= (-1)^ic;
        
        %SIGMA31 EN 0
        A_DWN((ic-2)*6+3+1,(ic-1)*6+1,:)=mic*S31(1,:);                     %Amplitude trans sens +
        A_DWN((ic-2)*6+3+1,(ic-1)*6+2,:)=mic*S31(2,:);                     %Amplitude longi sens +
        A_DWN((ic-2)*6+3+1,(ic-1)*6+3,:)=mic*S31(3,:);                     %Amplitude trans sens +
        
        %SIGMA32 EN 0
        A_DWN((ic-2)*6+3+2,(ic-1)*6+1,:)=mic*S32(1,:);                     %Amplitude trans sens +
        A_DWN((ic-2)*6+3+2,(ic-1)*6+2,:)=mic*S32(2,:);                     %Amplitude longi sens +
        A_DWN((ic-2)*6+3+2,(ic-1)*6+3,:)=mic*S32(3,:);                     %Amplitude trans sens +
        
        %SIGMA33 EN 0
        A_DWN((ic-2)*6+3+3,(ic-1)*6+1,:)=mic*S33(1,:);                     %Amplitude trans sens +
        A_DWN((ic-2)*6+3+3,(ic-1)*6+2,:)=mic*S33(2,:);                     %Amplitude longi sens +
        A_DWN((ic-2)*6+3+3,(ic-1)*6+3,:)=mic*S33(3,:);                     %Amplitude trans sens +
        
        %Amplitude des ondes avec kz sens +
        %U1 EN 0(ic)
        A_DWN((ic-2)*6+3+4,(ic-1)*6+1,:)= mic*u1(1,:,ic);                  %Amplitude trans sens +
        A_DWN((ic-2)*6+3+4,(ic-1)*6+2,:)= mic*u1(2,:,ic);                  %Amplitude longi sens +
        A_DWN((ic-2)*6+3+4,(ic-1)*6+3,:)= mic*u1(3,:,ic);                  %Amplitude trans sens +
        
        %U2 EN 0(ic)
        A_DWN((ic-2)*6+3+5,(ic-1)*6+1,:)= mic*u2(1,:,ic);                  %Amplitude trans sens +
        A_DWN((ic-2)*6+3+5,(ic-1)*6+2,:)= mic*u2(2,:,ic);                  %Amplitude longi sens +
        A_DWN((ic-2)*6+3+5,(ic-1)*6+3,:)= mic*u2(3,:,ic);                  %Amplitude trans sens +
        
        %U3 EN 0(ic)
        A_DWN((ic-2)*6+3+6,(ic-1)*6+1,:)= mic*u3(1,:,ic);                  %Amplitude trans sens +
        A_DWN((ic-2)*6+3+6,(ic-1)*6+2,:)= mic*u3(2,:,ic);                  %Amplitude longi sens +
        A_DWN((ic-2)*6+3+6,(ic-1)*6+3,:)= mic*u3(3,:,ic);                  %Amplitude trans sens +
        

    end
end

DWN.A_DWN = A_DWN;