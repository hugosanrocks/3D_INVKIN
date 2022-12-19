function [Gpr,Gsr]=Gij_P_SV_R(ks,kp,x,z,rij,g,C,dr)

% h0P = besselh(0,2,kp*rij);
% h0S = besselh(0,2,ks*rij);
% h2P = besselh(2,2,kp*rij);
% h2S = besselh(2,2,ks*rij);
n   = length(rij);
h0Pr= zeros(1,n);
h0Sr= zeros(1,n);
h2Pr= zeros(1,n);
h2Sr= zeros(1,n);

for k=1:n
    if x(k)==0 || abs(abs(x(k))-abs(z(k)))<1e-3
        h0Pr(k)=besselh_R_Im2(0,2,kp,0,rij(k),dr(k));
        h0Sr(k)=besselh_R_Im2(0,2,ks,0,rij(k),dr(k));
        h2Pr(k)=besselh_R_Im2(2,2,kp,0,rij(k),dr(k));
        h2Sr(k)=besselh_R_Im2(2,2,ks,0,rij(k),dr(k));
%     elseif abs(z(k))<1e-2*abs(x(k))
%         z1=-1e-2*abs(x(k))*sign(z(k));
%         h0Pr(k)=besselh_R_Im3(0,2,kp,x(k),z1,dr(k));
%         h0Sr(k)=besselh_R_Im3(0,2,ks,x(k),z1,dr(k));
%         h2Pr(k)=besselh_R_Im3(2,2,kp,x(k),z1,dr(k));
%         h2Sr(k)=besselh_R_Im3(2,2,ks,x(k),z1,dr(k));
    else
        h0Pr(k)=besselh_R_Im3(0,2,kp,x(k),z(k),dr(k));
        h0Sr(k)=besselh_R_Im3(0,2,ks,x(k),z(k),dr(k));
        h2Pr(k)=besselh_R_Im3(2,2,kp,x(k),z(k),dr(k));
        h2Sr(k)=besselh_R_Im3(2,2,ks,x(k),z(k),dr(k));
    end
end


% A   = h0P/C(1,1)+h0S/C(6,6);
% B   = h2P/C(1,1)-h2S/C(6,6);

% G   = zeros(2,2,n);
d   = eye(2);
%
% i=1;
% for j=1:2
%     G(i,j,:)=1/(1i*8)*(d(i,j)*A-(2*g(i,:).*g(j,:)-d(i,j)).*B);
% end
% 
% i=2;j=2;
% G(i,j,:)=1/(1i*8)*(d(i,j)*A-(2*g(i,:).*g(j,:)-d(i,j)).*B);
% G(2,1,:)=G(1,2,:);
% 
% end



Gpr	= zeros(2,2,n);
Gsr	= zeros(2,2,n);
% Gi	= zeros(2,2,n);

%longitudinal part
Ap  = h0Pr/C(1,1);
Bp  = h2Pr/C(1,1);

i=1;
for j=1:2
    Gpr(i,j,:)=1/(1i*8)*(d(i,j)*Ap-(2*g(i,:).*g(j,:)-d(i,j)).*Bp);
end
i=2;j=2;
Gpr(i,j,:)=1/(1i*8)*(d(i,j)*Ap-(2*g(i,:).*g(j,:)-d(i,j)).*Bp);
Gpr(2,1,:)=Gpr(1,2,:);

%shear part
As  = h0Sr/C(6,6);
Bs  =-h2Sr/C(6,6);

i=1;
for j=1:2
    Gsr(i,j,:)=1/(1i*8)*(d(i,j)*As-(2*g(i,:).*g(j,:)-d(i,j)).*Bs);
end
i=2;j=2;
Gsr(i,j,:)=1/(1i*8)*(d(i,j)*As-(2*g(i,:).*g(j,:)-d(i,j)).*Bs);
Gsr(2,1,:)=Gsr(1,2,:);

% %evanescente part
% Ai  = h0Pi/C(1,1)+h0Si/C(6,6);
% Bi  = h2Pi/C(1,1)-h2Si/C(6,6);
% 
% i=1;
% for j=1:2
%     Gi(i,j,:)=1/(1i*8)*(d(i,j)*Ai-(2*g(i,:).*g(j,:)-d(i,j)).*Bi);
% end
% i=2;j=2;
% Gi(i,j,:)=1/(1i*8)*(d(i,j)*Ai-(2*g(i,:).*g(j,:)-d(i,j)).*Bi);
% Gi(2,1,:)=Gi(1,2,:);

% G=Gpr+Gsr+Gi;
