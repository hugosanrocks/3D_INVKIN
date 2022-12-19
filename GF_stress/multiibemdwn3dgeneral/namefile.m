function name=namefile(para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creacion del nombre por plataforma %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strdim          ={'2D','2.5D','axi3D','gen3D'};
strpol          ={'SH','PSV'};
strfuente       ={'OP','FP'};
strmetodomed1   ={'FS','HS','DWN'};

if para.dim==1
    name=[strdim{para.dim},'_',strpol{para.pol},'_',strfuente{para.fuente}, ...
        '_xs',num2str(para.xs(1)),'_zs',num2str(para.zs(1)),'_npplo',num2str(para.npplo),'_'];
    if para.geo(1)==3 && para.nsubmed==1 && para.pol==1
        name=[name,'imagen','_'];
    else
        name=[name,strmetodomed1{para.geo(1)},'_'];
    end
else
    name=[strdim{para.dim},'_',strfuente{para.fuente}, ...
        '_xs',num2str(para.xs(1)),'_zs',num2str(para.zs(1)),'_npplo',num2str(para.npplo),'_'];
    if para.geo(1)==3 
        name=[name,strmetodomed1{para.geo(1)},'_'];
    end
end

if isunix
    name=[para.nomrep,'/',name,datestr(now,30)];
else
    name=[para.nomrep,'\',name,datestr(now,30)];
end
name=[name,'.mat'];