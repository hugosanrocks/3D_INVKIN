function [S1,S2]=S_PSV(rij,g,ks,kp,C,fij)
% trata el calculo de la contribution de un segmento centrado en xj,zj
% sobre un punto de colocacion xr,zr cual esta mas cercano que dr a xj,zj

% tij(I,J,:) : attention les i,j sont differents
% les minuscules i,j sont pour l IBEM : la contribution du segment j sur le
% point de colocation i
% les majuscules I,J, sont pour les tractions dans la direction I imposee
% par une force oriente selon J
% les tractions sont donc calculees tes que tnI=t(I,J).Phi_j
% et les tractions t(I,J) = S(I,1)*vn(1)+S(I,2)*vn(2)
% tn1 = squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:);
% tn2 = squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:);
% Aijx=[squeeze(tij(1,1,:));squeeze(tij(1,2,:))];%t_x=t_xj.Phi_j
% Aijz=[squeeze(tij(2,1,:));squeeze(tij(2,2,:))];%t_z=t_zj.Phi_j

% Ref: Sánchez-Sesma & Campillo 1991
% S1,S2 stress give a force in directio 1 or 2 resp.
% S#(1,1,:) sxx
% S#(1,2,:) sxz
% S#(2,2,:) szz

%version MP
S1  = zeros(2,2,length(rij));
S2  = zeros(2,2,length(rij));

h2P = besselh(2,2,kp*rij);
h2S = besselh(2,2,ks*rij);

B   = h2P/C(1,1)-h2S/C(6,6);

DP  = kp*rij.*besselh(1,2,kp*rij);
DS  = ks*rij.*besselh(1,2,ks*rij);

C0  = DP/C(1,1)-DS/C(6,6);

F1  = B+C(1,2)/(2*C(6,6)*C(1,1))*DP;
F2  = B+1/(2*C(6,6))*DS;
F3  = C0-4*B;

fac = 1i*C(6,6)./(2*rij);
if fij(1)~=0
    %f1
    %T(1,1) : force appliquee selon 1 T1=squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:);
    S1(1,1,:)=fac.*g(1,:).*( F1 + F2.*2 + F3.* g(1,:).^2);%.* vn(1,:);
    S1(1,2,:)=fac.*g(2,:).*( F2 + F3.* g(1,:).^2);%.*vn(2,:);
    
    %T(2,1) : force appliquee selon 1 T2=squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:);
    % S(1,2,:)=fac.*g(2,:).*( F2 + F3.*g(1,:).^2);%.*vn(1,:);
    S1(2,2,:)=fac.*g(1,:).*( F1 + F3.* g(2,:).^2);%.*vn(2,:);
end
if fij(2)~=0
    % T(1,2,:) : force appliquee selon 2 T1=squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:);
    S2(1,1,:)=fac.*g(2,:).*( F1 + F3.* g(1,:).^2);%.*vn(1,:);
    S2(1,2,:)=fac.*g(1,:).*( F2 + F3.* g(2,:).^2);%.*vn(2,:);
    
    % T(2,2,:) : force appliquee selon 2 T2=squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:);
    % S(1,2,:)=fac.*g(1,:).*( F2 + F3.* g(2,:).^2);%.*vn(1,:);
    S2(2,2,:)=fac.*g(2,:).*( F1 + 2*F2 + F3.* g(2,:).^2);%.*vn(2,:);
end