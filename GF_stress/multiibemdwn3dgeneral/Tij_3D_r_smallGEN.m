function [Tij,Sijf]=Tij_3D_r_smallGEN(coord,ii,jjpb,ks,kp,para,C,xr)
% Función de Green de TRACCIONES en el campo cercano
% trata el calculo de la contribution de un segmento centrado en xj,yj,zj
% sobre un punto de colocacion xi,yi,zi cual esta mas cercano que dr a xj,yj,zj
% La cubatura se hace sobre un elemento triangular.
%
% coord   :  coordenas de putnos de colucación / y fuentes virtuales
% ii      :  indice de l element (de surface) sur lequel se calcul la contribution de toutes
% jjpb    :  les autres sources virtuelles j (ponctuelle)
% ksi
% kpi

% #ij     :  Las coordenadas de los puntos de integración ya fueron calculadas en
%            la función malla_geom3Dgen.m con la función XYZ_FromBarycentric_Triang.m

nj      = length(jjpb);   % cantidad de fuentes virtuales
ngau    = para.gaussian.ngau;% cantidad de puntos gaussianos por fuente
Tij=zeros(3,3,nj);        % El resultado
Sijf=zeros(3,3,3,nj);
% El receptor:
% TODO : en lugar de ii, recibir los valores xr,yr,zr,vn
if nargin == 7
xr(1,1) = coord.x(ii);
xr(2,1) = coord.y(ii);
xr(3,1) = coord.z(ii);
vn(1,1) = coord.vnx(ii);
vn(2,1) = coord.vny(ii);
vn(3,1) = coord.vnz(ii);
elseif nargin == 8
% disp('Tij_3D_r_smallGEN: l30');disp(xr)
vn(1,1) = 0;
vn(2,1) = 0;
vn(3,1) = 1;
else
  error('faltan argumentos')
end
% Para cada fuente virtual se integre de todos los puntos
for iXi = 1:size(jjpb,2)
  % coordenadas y pesos de los putnos de la cubatura
  CubPts = coord.CubPts(1:3,1:ngau,jjpb(iXi));
  CubWts = para.cubature(1:ngau,3);
  
  % cosenos directores entre el receptor y las fuentes
  xij = xr(1)-CubPts(1,:);
  yij = xr(2)-CubPts(2,:);
  zij = xr(3)-CubPts(3,:);
  rij     = sqrt(xij.^2+yij.^2+zij.^2);
  g(1,:)  = xij./rij;
  g(2,:)  = yij./rij;
  g(3,:)  = zij./rij;
  
  %% Función de Green fuente puntual:
  [Trij,~,srikj] = Tij_3D(ks,kp,rij,g,C,vn);  %0,0
  indpb   = (rij==0);
  Trij(:,:,indpb) =0;
  srikj(:,:,:,indpb) =0;
  % suma de contribuciones de todos los puntos Gaussianos
  if nj~=1
    for i=1:3
      for j=1:3
        Tij(i,j,iXi)	= sum(CubWts.*squeeze(Trij(i,j,:)),1);
      end
    end
  else
    for i=1:3
      for j=1:3
        Tij(i,j,iXi)	= sum(CubWts.*squeeze(Trij(i,j,:)));
        for k=1:3
          Sijf(i,j,k,iXi) = sum(CubWts.*squeeze(srikj(i,j,k,:)));
        end
      end
    end
  end
end
% disp('here')
% return
end
