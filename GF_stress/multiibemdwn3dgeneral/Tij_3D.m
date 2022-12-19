function [T,TR,sikj]=Tij_3D(ks,kp,rij,gam,~,vnx,v_i,v_R)
% Función de Green de tracciones. Se usa teorema de reciprocidad.
% BSSA, Vol. 85, No. 1, pp. 269-284, 1995
% Seismic Response of Three-Dimensional Alluvial Valleys for Incident P, S, and Rayleigh Waves
% by Francisco J. Sanchez-Sesma and Francisco Luzon

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


T   = zeros(3,3,n);
d   = eye(3);
g0  = g(1,:)-g(2,:)-2*g(3,:);
fac = 1./(4*pi*rij.^2);

if nargin == 6
vn = vnx;
gknk= gam(1,:).*vn(1,:)+gam(2,:).*vn(2,:)+gam(3,:).*vn(3,:);
% Tracciones: (cálculo directo)
for i=1:3
    for j=1:3
        T(i,j,:)=fac.*(g0.*gam(i,:).*gam(j,:).*gknk +...
        g(3,:).*(gam(i,:).*vn(j,:)+gknk.*d(i,j))+g(2,:).*gam(j,:).*vn(i,:));
    end
end
TR = 0; sikj=0;
else
  % Esfuerzos: (el cálculo se puede reusar con el teorema de reciprocidad
  nR   = sum(v_R);
  TR  = zeros(3,3,nR);
  sikj = zeros(3,3,3,n); % tensor de esfuerzo, dada cada una de las fuentes
  for j = 1:3 % para cada dirección de la fuerza aplicada
    for i = 1:3
      for k = 1:3
        sikj(i,k,j,:) = fac.*(g0.*gam(i,:).*gam(j,:).*gam(k,:) + ...
          g(3,:).*(gam(i,:).*d(k,j)+gam(k,:).*d(i,j))+g(2,:).*gam(j,:).*d(k,i));
      end
      if v_R==0; continue; end
%     end
%     
%     for i = 1:3
    % Tracciones en el receptor (x) por todas las fuentes (xi)
    T(i,j,:) =  sikj(i,1,j,:)*vnx(1,v_i)...
               +sikj(i,2,j,:)*vnx(2,v_i)...
               +sikj(i,3,j,:)*vnx(3,v_i);
%   T(2,j,:) =  sikj(2,1,j,:)*vn(1)+sikj(2,2,j,:)*vn(2)+sikj(2,3,j,:)*vn(3);
%   T(3,j,:) =  sikj(3,1,j,:)*vn(1)+sikj(3,2,j,:)*vn(2)+sikj(3,3,j,:)*vn(3);
    
    % Tracciones en los receptores (xi) por la fuente en (x)
    TR(i,j,:) =  -squeeze(sikj(i,1,j,:)).*(vnx(1,v_R).')...
                 -squeeze(sikj(i,2,j,:)).*(vnx(2,v_R).')...
                 -squeeze(sikj(i,3,j,:)).*(vnx(3,v_R).');
    end
  end
end
% disp(max(max(max(T/TR)))) 
end