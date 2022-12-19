function Gij=Gij_3D_r_smallGEN(coord,xr,jjpb,ks,kp,para,C,ex)
% Función de Green de DESPLAZAMIENTO en el campo cercano
% trata el calculo de la contribution de un segmento centrado en xj,yj,zj
% sobre un punto de colocacion xi,yi,zi cual esta mas cercano que dr a xj,yj,zj
% La cubatura se hace sobre un elemento triangular.
%
% coord   :  coordenas de putnos de colucación / y fuentes virtuales
% xr yr zr:  punto o element (de surface) sur lequel se calcul la contribution de toutes
% jjpb    :  les autres sources virtuelles j (ponctuelle)
% ksi
% kpi

% #ij     :  Las coordenadas de los puntos de integración ya fueron calculadas en
%            la función malla_geom3Dgen.m con la función XYZ_FromBarycentric_Triang.m

nj      = length(jjpb);   % cantidad de fuentes virtuales
Gij     = zeros(3,3,nj);        % El resultado

% El receptor ya está dado en xr, yr, zr

% Para cada fuente virtual se integre de todos los puntos
for iXi = 1:size(jjpb,2)
  % coordenadas y pesos de los putnos de la cubatura
%   if ex
%   ngau    = para.gaussian.ngauex;% cantidad de puntos gaussianos por fuente
%   CubPts = coord.CubPtsex(1:3,1:ngau,jjpb(iXi));
%   CubWts = para.cubatureex(1:ngau,3);
%   else
  ngau    = para.gaussian.ngau;% cantidad de puntos gaussianos por fuente
  CubPts = coord.CubPts(1:3,1:ngau,jjpb(iXi));
  CubWts = para.cubature(1:ngau,3);
%   end
  % cosenos directores entre el receptor y las fuentes
  xij = xr(1)-CubPts(1,:);
  yij = xr(2)-CubPts(2,:);
  zij = xr(3)-CubPts(3,:);
  rij     = sqrt(xij.^2+yij.^2+zij.^2);
  g(1,:)  = xij./rij;
  g(2,:)  = yij./rij;
  g(3,:)  = zij./rij;
  
  %% Función de Green fuente puntual:
  grij    = Gij_3D(ks,kp,rij,g,C,length(rij));
  indpb   = (rij<=1E-5);
  grij(:,:,indpb) =0;
  % suma de contribuciones de todos los puntos Gaussianos
  if nj~=1
    for i=1:3
      for j=1:3
        Gij(i,j,iXi)	= sum(CubWts.*squeeze(grij(i,j,:)),1);
        Gij(j,i,iXi)	= Gij(i,j,iXi);
      end
    end
  else
    for i=1:3
      for j=1:3
        Gij(i,j,iXi)	= sum(CubWts.*squeeze(grij(i,j,:)));
        Gij(j,i,iXi)	= Gij(i,j,iXi);
      end
    end
  end
end
end