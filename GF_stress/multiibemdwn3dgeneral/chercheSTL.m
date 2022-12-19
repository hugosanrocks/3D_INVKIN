function [ nomrep ] = chercheSTL(imed)
[FileName,PathName] = uigetfile({'*.STL;*.stl;*.txt','ASCII Standard Triangle Language'}...
                      ,['Abrir archivo STL con geometría del medio ',num2str(imed)]...
                      ,'..'); % el el directorio padre
if isequal(FileName,0)
   disp('User selected Cancel')
   nomrep = '';
else
   nomrep = fullfile(PathName, FileName);
   disp(['User selected ', nomrep])
end
end

