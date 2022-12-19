% recargar listas
%%
for m = 2:para.nmed % para cada medio (sólo inclusiones)
  for p = 1:para.cont(m,1).NumPieces % cada pieza
    if (size(para.cont(m,1).piece{p}.fileName,2)>1) % se cargó algo válido
      if size(para.cont(m,1).piece{p}.geoFileData,2)>0 % y hay datos
        kind = para.cont(m,1).piece{p}.kind;
        if kind ~= 3 % si no es una frontera auxiliar
          if para.cont(m,1).piece{p}.isalist
            
            para.cont(m,1).piece{p} = loadStageList(...
                para.cont(m,1).piece{p},para.nf,para.cont(m,1).piece{p}.fileName);
            
          end
        end
      end
    end
  end
end