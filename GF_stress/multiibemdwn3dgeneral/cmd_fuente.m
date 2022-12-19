if para.fuente==1 %OP
    set(info.xzs        ,'string' ,'Origen fases OP');
    % tipo de onda
    if para.dim==1
        if para.pol==1
            if para.geo(1)==3
                set(info.tipo_onda,'visible','on');
                strtipo_onda ={'SH','L'};
                set(bouton.tipo_onda,'string',strtipo_onda,'visible','on');
            else
                set(info.tipo_onda,'visible','off');
                set(bouton.tipo_onda,'visible','off');
            end
        else
            set(info.tipo_onda,'visible','on');
            strtipo_onda ={'P','SV','R'};
            set(bouton.tipo_onda,'visible','on','string',strtipo_onda);
        end
    elseif para.dim>=3
        set(info.tipo_onda,'visible','on');
        strtipo_onda ={'P','SV','SH','R','L'};
        set(bouton.tipo_onda,'visible','on','string',strtipo_onda);
    end    
    set(info.gam        ,'visible','on');
    set(bouton.gam      ,'visible','on');
    set(info.orient     ,'visible','on');
else%FP
    set(info.xzs        ,'string' ,'Pos. Fuente');
    
    set(  info.tipo_onda,'visible','off');
    set(bouton.tipo_onda,'visible','off');
    
    if para.pol==1 && para.dim==1
        set(info.gam        ,'visible','off');
        set(bouton.gam      ,'visible','off');
        set(info.orient     ,'visible','off');
    else
        set(info.gam        ,'visible','on');
        set(bouton.gam      ,'visible','on');
        set(info.orient     ,'visible','on');
    end
end