function para = malla_geom3Dgen(para,fj)
% Refinar superficies para cumplir la número de puntos por longitud de onda
% Asume que la frecuencia va aumentando.
% Refina sólo los triángulos cuya área es tal que un círculo de la misma área
% y rádio 'r' es mayor que la 1/2 lambda/npplo. Donde lambda es la longitud
% de onda más pequeña en esa frontera.
%
% Subdivisión: Los triángulos se parten en dos desde la mitad del lado más
% largo al vértice opuesto. Se modi

%ya estamos en el ciclo de frecuencias
fj      = abs(fj);
npplo   = para.npplo;
%%%%%%%%%%%%%%%%%
%   subdividir  %
%%%%%%%%%%%%%%%%%
 maxSteps = 20; 

for m = 2:para.nmed % para cada medio (sólo inclusiones)
  for p = 1:para.cont(m,1).NumPieces % cada pieza
    kind = para.cont(m,1).piece{p}.kind;
    if kind ~= 3 % si no es una frontera auxiliar
      lambda  = para.cont(m,1).piece{p}.subdibData.minVel/fj; % longitud de onda más pequeña
      lambdaObjetivo = lambda/npplo; %('radio' máximo objetivo)
      radioObjetivo = lambdaObjetivo / 2;
%               disp(['lambda Objetivo ',num2str(lambdaObjetivo)])
%               disp(['radio  Objetivo ',num2str(radioObjetivo)])
%               disp(['maxRadio ',num2str(max(para.cont(m,1).piece{p}.subdibData.radios))])
      % indices de los triángulos que se deben subdividir:
      deNuevo = true; iStep = 0; 
      while deNuevo
%         l_faltantes = []; 
        iStep = iStep + 1;
        l = find(para.cont(m,1).piece{p}.subdibData.radios>radioObjetivo);
        l = reshape(l,1,length(l));
        %         disp(['at ' num2str(fj) ' on p',num2str(p),' splitting ' num2str(length(l))]);
        if ~isempty(l)
          NewTriangs = [];
          for t = l % cada triángulo grandote se parte en dos
            twoNewTriangs = splitThisTriangle(para.cont(m,1).piece{p}.subdibData,t);
            NewTriangs = apilar(NewTriangs,twoNewTriangs);
          end
          % borrar los triángulos l
          para.cont(m,1).piece{p}.subdibData = borrar(l,para.cont(m,1).piece{p}.subdibData);
          
          % agregar los nuevos triángulos
          para.cont(m,1).piece{p}.subdibData = apilar(NewTriangs,para.cont(m,1).piece{p}.subdibData);
          
          disp(['m',num2str(m),' p',num2str(p),' newT_',num2str(length(NewTriangs.centers))])
        else
          deNuevo = false;
        end
        if iStep > maxSteps; warning('max subdiv steps reached'); deNuevo = false; end
      end % subdividir deNuevo
%       disp(['m',num2str(m),' p',num2str(p),' nColocPts:_', ...
%         num2str( size(para.cont(m,1).piece{p}.subdibData.centers,2))]); %cant. triangs
    end
  end
end


% Puntos de cubatura
for m = 2:para.nmed % para cada medio (sólo inclusiones)
  for p = 1:para.cont(m,1).NumPieces % cada pieza
    kind = para.cont(m,1).piece{p}.kind;
    if kind ~= 3 % si no es una frontera auxiliar
      [para.cont(m,1).piece{p}.subdibData.trianglesCubaturePoints]...
        = XYZ_FromBarycentric_Triang(...
        para.cont(m,1).piece{p}.subdibData.triangles,...
        para.cubature);
%       [para.cont(m,1).piece{p}.subdibData.trianglesCubaturePointsex]...
%         = XYZ_FromBarycentric_Triang(...
%         para.cont(m,1).piece{p}.subdibData.triangles,...
%         para.cubatureex);
      
      %         %% Test me:
      %         figure(343); hold on
      %         sdD = para.cont(m,1).piece{p}.subdibData;
      %         for t = 1:size(sdD.triangles,3)
      %           plot3(sdD.triangles(1,[1 2 3 1],t),...
      %             sdD.triangles(2,[1 2 3 1],t),...
      %             sdD.triangles(3,[1 2 3 1],t),'k-','LineWidth',0.5)
      %           plot3(sdD.centers(1,t),...
      %                 sdD.centers(2,t),...
      %                 sdD.centers(3,t),'k.')
      %           plot3(sdD.trianglesCubaturePoints(1,:,t),...
      %                 sdD.trianglesCubaturePoints(2,:,t),...
      %                 sdD.trianglesCubaturePoints(3,:,t),'r.')
      %         end
      %         clear sdD
      %         %%
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables para ensamblar dominios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbpt0 = 0; % cantidad de puntos de colocación
medio(para.nmed).ind =[];
for m=2:para.nmed;
  for p = 1:para.cont(m,1).NumPieces
    kind = para.cont(m,1).piece{p}.kind;
    if kind ~= 3 % si no es una frontera auxiliar
      nbpt0 = nbpt0 + size(para.cont(m,1).piece{p}.subdibData.centers,2);
    end
  end
end
%% inicio
n1          = para.nmed;
coord       = struct( ...
  'x'     ,zeros(1,nbpt0),	'y'     ,zeros(1,nbpt0) , 'z'   ,zeros(1,nbpt0), ...
  'vnx'   ,zeros(1,nbpt0),	'vny'   ,zeros(1,nbpt0) , 'vnz' ,zeros(1,nbpt0), ...
  'drxz'  ,zeros(1,nbpt0),  'cv'    ,zeros(nbpt0,1) , 'dA'  ,zeros(1,nbpt0),...
  'Xim' 	,zeros(nbpt0,n1), 'phi'  	,zeros(nbpt0,n1), ...
  'CubPts',zeros(3,para.gaussian.ngau,nbpt0)           , ...
  'fl'  	,zeros(nbpt0,n1), 'indm'  ,[]             , ...
  'ju', [], 'nbeq', [], 'nfl', [],...
  'm',zeros(1,nbpt0),'p',zeros(1,nbpt0),...
  'mint',zeros(1,nbpt0),'mext',zeros(1,nbpt0));%, ...
%   'x0',zeros(1,nbpt0) ,'th',zeros(1,nbpt0), 'nbptc',zeros(1,nbpt0));
coord.nbpt  = nbpt0;
%% x,y,z,vnx,vny,vnz,cv,Xim,phi,fl;  (no r)
jjF = 0; jp  = 0;    %initializacion del indice de los phi
for m=2:para.nmed;
  for p = 1:para.cont(m,1).NumPieces
    kind = para.cont(m,1).piece{p}.kind;
    if kind ~= 3 % si no es una frontera auxiliar
      n = size(para.cont(m,1).piece{p}.subdibData.centers,2); %cant. triangs
      jjI = jjF + 1;
      jjF = jjF + n;
      jj = jjI:jjF; %indices de putnos de colocación actuales
      
      %indices para recordar de donde salió
      coord.m(jj) = m;
      coord.p(jj) = p;
      % lo demás:
      coord.x(jj) = para.cont(m,1).piece{p}.subdibData.centers(1,:);
      coord.y(jj) = para.cont(m,1).piece{p}.subdibData.centers(2,:);
      coord.z(jj) = para.cont(m,1).piece{p}.subdibData.centers(3,:);
      %     coord.x0(jj) =
      coord.vnx(jj) = para.cont(m,1).piece{p}.subdibData.N(:,1);
      coord.vny(jj) = para.cont(m,1).piece{p}.subdibData.N(:,2);
      coord.vnz(jj) = para.cont(m,1).piece{p}.subdibData.N(:,3);
      coord.cv(jj) = para.cont(m,1).piece{p}.subdibData.cv;
      
      % puntos de integración Gaussiana
      coord.CubPts  (:,:,jj) = para.cont(m,1).piece{p}.subdibData.trianglesCubaturePoints  (:,:,:);
%       coord.CubPtsex(:,:,jj) = para.cont(m,1).piece{p}.subdibData.trianglesCubaturePointsex(:,:,:);
      % de la inclusión:
        mext = para.cont(m,1).piece{p}.continuosTo; 
        coord.mext(jj) = mext; %medio exterior
        mint = m;
        coord.mint(jj) = mint;%medio interior
%       mext	= para.cont1(m-1).vec.mv; coord.mext(jj) = mext; %medio exterior
%       mint	= para.cont1(m-1).m;      coord.mint(jj) = mint;%medio interior
      c     = coord.cv(jj);     %1 contorno de arriba; 2 contorno de abajo
      
      %indexacion (ver si hay un punto o dos), medios en contacto, y normal
      if mext==0 || para.reg(mext).rho==0 %2015-01-16
        %superficie libre y hueco afuera
        coord.Xim(jj,mint)  = c;	%ademas de saber que el punto Xi pertenece al medio m, se conoce su posicion arriba/abajo en el contorno
        coord.phi(jj,mint)  = jp+(1:n); %posicion de los phi en el vector X del sistema: AX=B
        jp                  = jp+n;
      elseif mint==0 || para.reg(m).rho==0
        %para tomar en cuenta los huecos
        % continuidad de tracciones:
        coord.Xim(jj,mext)  = (c==2)+2*(c==1);%inversion de la posicion
        coord.phi(jj,mext)  = jp+(1:n);
        jp                  = jp+n;
      else %continuidad
        % de tracciones:
        coord.Xim(jj,mext)  = (c==2)+2*(c==1);%inversion de la posicion
        coord.phi(jj,mext)  = jp+(1:n);
        jp                  = jp+n;
        if kind == 2 % si es una frontera de continuidad
          % de desplazamientos:
          coord.Xim(jj,mint)	= c;
          coord.phi(jj,mint) 	= jp+(1:n);
          jp                  = jp+n;
        else
          mext=0;
          para.cont1(m-1).vec.mv = 0;
          coord.mext(jj) = 0;
        end
      end
      
      if mint~=0
        medio(mint).ind	= [medio(mint).ind;[jj(1) jj(end)]];%principio y fin
        if para.reg(mint).bet==0 && para.reg(mint).rho~=0
          %indicacion de presencia de fluido
          coord.fl(jj,mint)	= 1;
        end
      end
      if mext~=0
        medio(mext).ind	= [medio(mext).ind;[jj(1) jj(end)]];%principio y fin
        if  para.reg(mext).bet==0 && para.reg(mext).rho~=0
          %indicacion de presencia de fluido
          coord.fl(jj,mext)    = 1;
        end
      end
      
      %se censa los puntos que pertenecen al medio m o mas bien al sub medio mi
      for mi=1:size(coord.Xim,2)
        coord.indm(mi).ind=(coord.Xim(:,mi)~=0);
      end
      
      %se censa los puntos que no tienen frontera libre de manera a aplicar la
      %continuidad de los desplazamientos
      coord.ju        = cumsum(sum(coord.Xim~=0,2)==2);
      coord.nbeq      = coord.ju(coord.nbpt)+coord.nbpt;
      coord.nfl       = sum(sum(coord.fl));
      
      % ".drxy" longitud del arco que describe la geometría axisimétrica
      % ".drxz" distancia entre centros de segmentos discretiación 2D
      
      % ".dA"   área del segmento 'escudo'.
      %         Se da el área del triángulo
      coord.dA(jj) = para.cont(m,1).piece{p}.subdibData.areas(:);
      
      %         Se substituye por el radio de un círculo de la misma área:
      coord.drxz(jj) = sqrt(coord.dA(jj)/pi);
      
    end
  end
end
para.coord = coord;
end

function out = apilar(a,b)
if isempty(a); out=b; return;end
out = a;
ni = size(a.centers,2)+1;
nf = ni-1 + size(b.centers,2);
l = ni:nf;
if ~isempty(l)
  out.triangles(:,:,l) = b.triangles;
  out.centers(:,l) = b.centers;
  out.areas(l) = b.areas;
  out.N(l,:) = b.N;
  out.cv(l) = b.cv;
  out.radios(l) = b.radios;
  %   out.minVel(l) = b.minVel;
end
end

function out = borrar(l,a)
out = a;
if ~isempty(l)
  out.triangles(:,:,l) = [];
  out.centers(:,l) = [];
  out.areas(l) = [];
  out.N(l,:) = [];
  out.cv(l) = [];
  out.radios(l) = [];
  % minVel   es un escalar válido para toda la pieza
end
end