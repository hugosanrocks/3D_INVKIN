function [medio] = inclusiontest3G(xs,ys,zs,para)
% El medio al que pertenece un punto(x,y,z) dado un polihedro en 3D
% Los contornos están formados por la rectificación de su superfice
% mediante triángulos. Los contornos son superfices cerradas

npts    = size(xs,1);
medio   = ones(npts,1); %Está en el medio de fondo (suposición inicial)
if para.geo(1) > 1 % Semiespacio homogeneo ó estratificado
  medio(zs<0) = 0; % Medio 0
end

m       = para.nmed;
dr=1e-6; % tolerancia
testme = ones(npts,1); % 1: si revisar

while m>1
  if isfield(para.cont(m,1),'FV')
    if ~isempty(para.cont(m,1).FV)
    for ip = 1:npts
     if xs(ip)-min(para.cont(m,1).FV.vertices(:,1))<-dr || max(para.cont(m,1).FV.vertices(:,1))-xs(ip)<-dr || ...
        ys(ip)-min(para.cont(m,1).FV.vertices(:,2))<-dr || max(para.cont(m,1).FV.vertices(:,2))-ys(ip)<-dr || ...
        zs(ip)-min(para.cont(m,1).FV.vertices(:,3))<-dr || max(para.cont(m,1).FV.vertices(:,3))-zs(ip)<-dr
      %de plano, no pertenece a m
      testme(ip)=0; % 0: no revisar
     end
    end
    t = find(testme); % indices de los que hay que revisar
    if ~isempty(t)
      in = inpolyhedron(para.cont(m,1).FV,...
                        [xs(t),ys(t),zs(t)]);%,'facenormals', para.cont(m,1).FV.facenormals(t,1:3));
      medio(t(in)) = m;
    end
    end
  end
    m=m-1;
end
end