if isfield(para,'gamr')
    para=rmfield(para,'gamr');
end
if isfield(para,'kxk')
    para=rmfield(para,'kxk');
end
if isfield(para,'kzsigno')
    para=rmfield(para,'kzsigno');
end
if isfield(para,'df')
    para=rmfield(para,'df');
end
if isfield(para,'cont0')
    para=rmfield(para,'cont0');
end
if isfield(para,'nmed0')
    para=rmfield(para,'nmed0');
end
if isfield(para,'subm1')
    para=rmfield(para,'subm1');
end
if isfield(para,'subm')
    para=rmfield(para,'subm');
end
if isfield(para,'fus')
    para=rmfield(para,'fus');
end
if isfield(para,'cont1')
    para=rmfield(para,'cont1');
end
if isfield(para,'nmed1')
    para=rmfield(para,'nmed1');
end
if isfield(para,'nmedf')
    para=rmfield(para,'nmedf');
end
if isfield(para,'listc')
    para=rmfield(para,'listc');
end
if isfield(para,'fuenteimagen')
    para=rmfield(para,'fuenteimagen');
end
if isfield(para,'xzs')
    para=rmfield(para,'xzs');
end


if isfield(para.cont,'long')
    para.cont=rmfield(para.cont,'long');
end
if isfield(para.cont,'vec')
    para.cont=rmfield(para.cont,'vec');
end
if isfield(para.cont,'int')
    para.cont=rmfield(para.cont,'int');
end
if isfield(para.cont,'seg')
    para.cont=rmfield(para.cont,'seg');
end
