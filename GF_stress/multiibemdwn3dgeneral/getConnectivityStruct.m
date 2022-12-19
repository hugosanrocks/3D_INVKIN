function para = getConnectivityStruct(para)
% Hacer struct de conectividades y geometría inicial
% Se devuelve cont1() con los elementos por cada contorno:
%     .m   :  El medio al que pertenece el contorno
%     .mv  :  El medio con el que hay continuidad de desp. y trac.
%     .vec :  Estructura con los siguientes elementos:
%         .xc  .yc  .zc   : componentes centro de elemento
%         .vnx .vny .vnz  : componentes normal orientada hacia -z
%                           Las normales se reorientan pero el orden de los
%                           vértices no se actualiza (aguas con crossprod!)
%         .cv             : [1 ó 2] si es contorno de arriba o abajo
% (no)    .r              : logitud rectificada acumulada (2D)
%         .radio          : círculo estimado del área calculada
%         .minVel         : la menor de las velocidades a cada lado del punto
n = 0;
for m = 2:para.nmed  % 1 ?
  n = n + para.cont(m,1).NumPieces;
end
cont1 = struct('m', zeros(1,n),...
              'mv', zeros(1,n),...
             'vec', [],...
           'radio', zeros(1,n),...
          'minVel', zeros(1,n));
i = 1;
hacerlo = false;
for m = 2:para.nmed
  for p = 1:para.cont(m,1).NumPieces
    if (size(para.cont(m,1).piece{p}.fileName,2)>1) % se cargó algo válido
      if size(para.cont(m,1).piece{p}.geoFileData,2)>0 % y hay datos
        kind = para.cont(m,1).piece{p}.kind;
        if kind ~= 3 % si no es una frontera auxiliar
            hacerlo = true;
          cont1(i).m = m; %medio al que pertenece
          cont1(i).vec.mv = para.cont(m,1).piece{p}.continuosTo;
          % por cada cara (triángulo)
          cont1(i).vec.xc = squeeze(para.cont(m,1).piece{p}.geoFileData.centers(1,:));
          cont1(i).vec.yc = squeeze(para.cont(m,1).piece{p}.geoFileData.centers(2,:));
          cont1(i).vec.zc = squeeze(para.cont(m,1).piece{p}.geoFileData.centers(3,:));
          cont1(i).vec.vnx = squeeze(para.cont(m,1).piece{p}.geoFileData.N(:,1));
          cont1(i).vec.vny = squeeze(para.cont(m,1).piece{p}.geoFileData.N(:,2));
          cont1(i).vec.vnz = squeeze(para.cont(m,1).piece{p}.geoFileData.N(:,3));
          range = cont1(i).vec.vnz > 0;
          cont1(i).vec.cv = ones(size(cont1(i).vec.xc)); % es contorno de arriba
          cont1(i).vec.cv(range) = 2; % es contorno de abajo
          % por convención en multi-ibem, las normales apuntan a -z
          % invertir la normal
          cont1(i).vec.vnx(range) = -cont1(i).vec.vnx(range);
          cont1(i).vec.vny(range) = -cont1(i).vec.vny(range);
          cont1(i).vec.vnz(range) = -cont1(i).vec.vnz(range);
          cont1(i).radio = sqrt(para.cont(m,1).piece{p}.geoFileData.areas(:,1)/pi);
          % la velocidad más baja de cada lado del contorno en el punto
          indm    = [cont1(i).m;cont1(i).vec.mv]; % los dos medios involucrados
          indm(indm==0)   = []; 
          indm    = squeeze(indm);
          v       = 0*indm; %init
          for k=1:length(indm)
            if para.reg(indm(k)).rho==0 % si está hueco
              v(k) = 0;
            else
                if isempty(para.reg(indm(k)).sub) %              si la regíon referida,
                   v(k) = para.reg(indm(k)).bet; %               es un medio cerrado o,
                else%                                            es fondo estratificado
                   v(k) = min(para.reg(indm(k)).sub(1:end).bet); 
                end
              if v(k)==0 % si es un medio acústico
                  if isempty(para.reg(indm(k)).sub)
                     v(k) = para.reg(indm(k)).alpha;
                  else
                     v(k) = min(para.reg(indm(k)).sub(1:end).alpha);
                  end
              end
            end
          end
          v(v==0) = [];
          v       = squeeze(v);
          if isempty(v)
            error('v isempty')
          end
          cont1(i).minVel = min(v);
          
          % En la estructura de datos original
          D = para.cont(m,1).piece{p}.geoFileData; 
          % se preservan los campos: .centers .triangles .areas .N
          D = rmfield(D,{'fileName','F','V'}); % estas no
          
          D.cv = ones(size(cont1(i).vec.xc)); % es contorno de arriba
          D.cv(range) = 2; % es contorno de abajo
          % por convención en multi-ibem, las normales apuntan a -z
          % invertir la normal
          D.N(range,1) = -D.N(range,1);
          D.N(range,2) = -D.N(range,2);
          D.N(range,3) = -D.N(range,3);
          D.radios = sqrt(D.areas(:,1)/pi);
          D.minVel = min(v);
          if isfield(para.cont(m,1).piece{p},'subdibData')
            para.cont(m,1).piece{p} = rmfield(para.cont(m,1).piece{p},'subdibData');
          end
          para.cont(m,1).piece{p}.subdibData = D;
          i = i + 1;
        end
      else
        warning('no todos las piezas tienen datos');
      end % hay datos
    end % se cargó
  end % cada pieza
end % cada medio
if  hacerlo
para.cont1 = cont1;

% de cont0 se toma la geometría para conocer la región de los receptores 
para.cont0  = cont1;
para.nmed0	= 1 + length(para.cont0);

%cuando hay contornos que se dividen, hay que reencontrar el medio exterior
%con la nueva numerotacion
nmed1  = length(cont1);
subm1(1) = 1; % el medio de fondo
for m1=1:nmed1 % las inclusiones
    if ~isempty(m1)
        subm1(m1+1)=cont1(m1).m;
    end
end
for m=1:para.nmed
    subm(m).m=find(subm1==m);
end
para.subm1=subm1;
para.subm =subm;
end %hacerlo
%% ...se va a indentificar los segmentos que son frontera entre 2 dos medios o de frontera libre
% j   = 1;


end