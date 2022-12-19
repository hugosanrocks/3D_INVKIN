set(bouton.ninc,'string',para.ninc);
if para.ninc~=length(para.gam)
    itmp    = get(bouton.inc,'value');
    xstmp   = str2double(get(bouton.xs ,'string'));
    zstmp   = str2double(get(bouton.zs ,'string'));

    para.strinc =[];
    para    = rmfield(para,'gam');
    para    = rmfield(para,'xs');
    para    = rmfield(para,'zs');
    para    = rmfield(para,'tipo_onda');
    if para.dim>=3
        ystmp   = str2double(get(bouton.ys ,'string'));
        para    = rmfield(para,'ys');
        para    = rmfield(para,'phi');
    end
    for i=1:para.ninc
        tmpstr=num2str(i);
        while length(tmpstr)<floor(log10(para.ninc)+1)
            tmpstr=['0',tmpstr];
        end
        para.strinc         = [para.strinc;tmpstr];
        para.xs(i)          = xstmp;
        para.zs(i)          = zstmp;
        para.tipo_onda(i)   = 1;
        if para.dim>=3
            para.ys(i)      = ystmp;
        end
    end
    para.gam  	= linspace(-90,90,para.ninc);
    para.phi  	= zeros(1,para.ninc);

    set(bouton.inc,'string',para.strinc,'value',1);
    if para.rafraichi==0
        cmd_inc;
    end
else
%     if para.rafraichi==0
%         dibujo_conf_geo(para,bouton.axe_conf_geo,bouton);
%     end
end

