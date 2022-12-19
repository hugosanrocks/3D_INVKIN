needsToUpdate = false;
for m = 2:para.nmed % para cada medio (sólo inclusiones)
  for p = 1:para.cont(m,1).NumPieces % cada pieza
    if (size(para.cont(m,1).piece{p}.fileName,2)>1) % se cargó algo válido
      if size(para.cont(m,1).piece{p}.geoFileData,2)>0 % y hay datos
        kind = para.cont(m,1).piece{p}.kind;
        if kind ~= 3 % si no es una frontera auxiliar
          if para.cont(m,1).piece{p}.isalist
            needsToUpdate = true;
            auxcont.fileName = para.cont(m,1).piece{p}.stage(end).fileName;
            if isfield(para.cont(m,1).piece{p},'geoFileData')
              para.cont(m,1).piece{p} = rmfield(para.cont(m,1).piece{p},'geoFileData');
            end
            if isfield(para.cont(m,1).piece{p},'subdibData')
              para.cont(m,1).piece{p} = rmfield(para.cont(m,1).piece{p},'subdibData');
            end
            [para.cont(m,1).piece{p}.geoFileData,flag] = ...
              previewSTL(0,auxcont);
            clear auxcont
          end
        end
      end
    end
  end
end
if(needsToUpdate)
  para = getConnectivityStruct(para);
end
clear m p kind needsToUpdate 