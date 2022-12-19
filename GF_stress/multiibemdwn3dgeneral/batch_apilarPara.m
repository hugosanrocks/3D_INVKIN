function batch_apilarPara(para)
% Copies all of 'para' into a new element of struct batch
% The intention is to run calculo(para) many times.
% If batch.mat does not exist it is created

  % borrar el archivo si es especificado
  para.name  = namefile(para);
cd ..
cd out
if ~isstruct(para)
  if para==false
    if exist(fullfile(cd, 'batch.mat'), 'file')
      disp('deleting batch.mat')
      delete('batch.mat')
    else
      disp('the file did not existed anyway')
    end
  end
else
if ~exist(fullfile(cd, 'batch.mat'), 'file')
  % si no existe el archivo lo creamos
  batch = cell(1);
  batch{1,1} = para;
  batch{1,2} = []; % result
  disp('saving fresh batch.mat')
  save('batch.mat','batch')
else
  % si existe, apilar
  load('batch.mat','batch') %cargar los originales,
  len = size(batch,1);
  batch{len+1,1} = para;
  batch{len+1,2} = []; % result
  disp('saving onto batch.mat')
  save('batch.mat','batch')
end
end
cd ..
cd multi-dwn-ibem.matlab
end


