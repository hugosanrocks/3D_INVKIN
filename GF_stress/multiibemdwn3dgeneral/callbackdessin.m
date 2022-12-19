function [RESULT,b_dib]=callbackdessin(para,nomcarpeta,bdessin,b_dib)

% if exist('vardesalida','var')==1
%     tmpvardesalida  = vardesalida;
% end

cd(para.nomrep)
[nomfile2,nomrep2]=uigetfile('*.mat','archivo de datos y salida');
if isequal(nomrep2,0) || isequal(nomfile2,0)
    errordlg('entra un archivo de tipo configparatmp.mat','info');
else
    name=[nomrep2,nomfile2];
    load(name);
    b_dib(1).name=name;
end
para.redraw=1;

cd(nomcarpeta)
if exist('utc','var')
    if ~exist('stc','var')
        stc = 0;
        sw  = 0;
    end
    if ~isfield(para,'spct')
        para.spct = 0;
    end
    b_dib=dibujo(para,bdessin,utc,uw,stc,sw,cont1,b_dib);
elseif exist('uw','var')
    utc=0;stc=0;
    b_dib=dibujo(para,bdessin,utc,uw,stc,sw,cont1,b_dib);
else
    utc=0;cont1=[];uw=0;stc=0;sw=0;name=[];
end

RESULT.utc=utc;
RESULT.uw=uw;
RESULT.stc=stc;
RESULT.sw=sw;
RESULT.name=name;
RESULT.cont1=con1;
end
