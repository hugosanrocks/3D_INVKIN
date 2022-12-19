function [S_fx,S_fy,S_fz]=S_3D(rij,gam,ks,kp,fij)
% trata el calculo de los esfuerzos, se basa sobre Tij_3D identificando las
% proyeciones

S_fx  = zeros(3,3,length(rij));
S_fy  = zeros(3,3,length(rij));
S_fz  = zeros(3,3,length(rij));


n       = length(rij);
ba      = kp/ks;%beta/alpha
kpr     = kp*rij;
ksr     = ks*rij;
ksrm1   = 1./ksr;
g       = zeros(3,n);

g(1,:)= ...
    (                      4-12i*ksrm1-12*ksrm1.^2).*exp(-1i*ksr)+...
    (-1i*ba*ksr-4*ba^2-1+ 12i*ba*ksrm1+12*ksrm1.^2).*exp(-1i*kpr);

g(2,:)= ...
    (                             -2+6i*ksrm1+6*ksrm1.^2).*exp(-1i*ksr)+...
    ( 1i*(2*ba^3-ba)*ksr+4*ba^2-1-6i*ba*ksrm1-6*ksrm1.^2).*exp(-1i*kpr);

g(3,:)= ...
    (-1i*ksr-3+ 6i*ksrm1+6*ksrm1.^2).*exp(-1i*ksr)+...
    (2*ba^2 -6i*ba*ksrm1-6*ksrm1.^2).*exp(-1i*kpr);


g0  = g(1,:)-g(2,:)-2*g(3,:);
fac = 1./(4*pi*rij.^2);

if fij(1)~=0
    %T(1,1) : force appliquee selon 1 T1=squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:)+squeeze(S(1,3,:)).'.*vn(3,:);
    S_fx(1,1,:) = fac.*(g0.*gam(1,:).^2 + 2*g(3,:)+g(2,:)).*gam(1,:);%.*vn(1,:);
    S_fx(1,2,:) = fac.*(g0.*gam(1,:).^2 +   g(3,:)       ).*gam(2,:);%.*vn(2,:);
    S_fx(1,3,:) = fac.*(g0.*gam(1,:).^2 +   g(3,:)       ).*gam(3,:);%.*vn(3,:);
    
    %T(2,1) : force appliquee selon 1 T2=squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:)+squeeze(S(3,2,:)).'.*vn(3,:);
    S_fx(2,1,:) = S_fx(1,2,:);
    S_fx(2,2,:) = fac.*(g0.*gam(2,:).^2 + g(2,:)).*gam(1,:);
    S_fx(2,3,:) = fac.* g0.*gam(1,:).*gam(2,:).*gam(3,:);
    
    %T(3,1) : force appliquee selon 1 T3=squeeze(S(1,3,:)).'.*vn(1,:)+squeeze(S(2,3,:)).'.*vn(2,:)+squeeze(S(3,3,:)).'.*vn(3,:);
    S_fx(3,1,:) = S_fx(1,3,:);
    S_fx(3,2,:) = S_fx(2,3,:);
    S_fx(3,3,:) = fac.*(g0.*gam(3,:).^2 + g(2,:)).*gam(1,:);
end
if fij(2)~=0
    %T(1,2) : force appliquee selon 2 T12=squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:)+squeeze(S(1,3,:)).'.*vn(3,:);
    S_fy(1,1,:) = fac.*(g0.*gam(1,:).^2 + g(2,:)).*gam(2,:);
    S_fy(1,2,:) = fac.*(g0.*gam(2,:).^2 + g(3,:)).*gam(1,:);
    S_fy(1,3,:) = fac.* g0.*gam(1,:).*gam(2,:)   .*gam(3,:);
    
    %T(2,2) : force appliquee selon 2 T22=squeeze(S(1,2,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:)+squeeze(S(3,2,:)).'.*vn(3,:);
    S_fy(2,1,:) = S_fy(1,2,:);
    S_fy(2,2,:) = fac.*(g0.*gam(2,:).^2 + g(2,:) + 2*g(3,:)).*gam(2,:);
    S_fy(2,3,:) = fac.*(g0.*gam(2,:).^2          +   g(3,:)).*gam(3,:);
    
    %T(3,2) : force appliquee selon 2 T32=squeeze(S(1,3,:)).'.*vn(1,:)+squeeze(S(2,3,:)).'.*vn(2,:)+squeeze(S(3,3,:)).'.*vn(3,:);
    S_fy(3,1,:) = S_fy(1,3,:);
    S_fy(3,2,:) = S_fy(2,3,:);
    S_fy(3,3,:) = fac.*(g0.*gam(3,:).^2 +g(2,:)).*gam(2,:);
    
    %verification des correspondances entre f_y et fx: y_y=x_x et x_y=-y_x
    %=> g(2)_fy=g(1)_fx et g(1)_fy=-g(2)_fx 
%     S_fy(1,1,:) = S_fx(2,2,:);
%     S_fy(1,2,:) =-S_fx(1,2,:);
%     S_fy(1,3,:) =-S_fx(2,3,:);
%     S_fy(2,1,:) = S_fy(1,2,:);
%     S_fy(2,2,:) = S_fx(1,1,:);
%     S_fy(2,3,:) = S_fx(1,3,:);
%     S_fy(3,1,:) = S_fy(1,3,:);
%     S_fy(3,2,:) = S_fy(2,3,:);
%     S_fy(3,3,:) = S_fx(3,3,:);
end
if fij(3)~=0
    %T(1,3) : force appliquee selon 3 T1=squeeze(S(1,1,:)).'.*vn(1,:)+squeeze(S(1,2,:)).'.*vn(2,:)+squeeze(S(1,3,:)).'.*vn(3,:);
    S_fz(1,1,:) = fac.*(g0.*gam(1,:).^2 + g(2,:)).*gam(3,:);
    S_fz(1,2,:) = fac.* g0.*gam(1,:) .*gam(2,:)  .*gam(3,:);
    S_fz(1,3,:) = fac.*(g0.*gam(3,:).^2 + g(3,:)).*gam(1,:);
    
    %T(2,3) : force appliquee selon 3 T2=squeeze(S(2,1,:)).'.*vn(1,:)+squeeze(S(2,2,:)).'.*vn(2,:)+squeeze(S(2,3,:)).'.*vn(3,:);
    S_fz(2,1,:) = S_fz(1,2,:);
    S_fz(2,2,:) = fac.*(g0.*gam(2,:).^2+g(2,:)).*gam(3,:);
    S_fz(2,3,:) = fac.*(g0.*gam(3,:).^2+g(3,:)).*gam(2,:);
    
    %T(3,3) : force appliquee selon 3 T3=squeeze(S(3,1,:)).'.*vn(1,:)+squeeze(S(3,2,:)).'.*vn(2,:)+squeeze(S(3,3,:)).'.*vn(3,:);
    S_fz(3,1,:) = S_fz(1,3,:);
    S_fz(3,2,:) = S_fz(2,3,:);
    S_fz(3,3,:) = fac.*(g0.*gam(3,:).^2 + g(2,:) + 2*g(3,:)).*gam(3,:);
end