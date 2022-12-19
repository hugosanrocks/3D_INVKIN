med     =get(bouton.med,'value');
if med==1 && para.geo(med) == 3
    return;
end
if para.chggeo==1
    para.chggeo=0;
    if para.geo(med) == 1 && med>1
        para.cont(med,1).ba   = 0.25;
        para.cont(med,2).ba   = 0.25;
    else
        para.cont(med,1).ba   = 0;
        para.cont(med,2).ba   = 0;
        if para.geo(1)==2
            para.cont(med,1).a =para.cont(1,1).a;
            para.cont(med,1).xa=para.cont(1,1).xa;
        else
            para.cont(med,1).a = 6;
            para.cont(med,1).xa=-3;
        end
    end
end
set(bouton.bet  ,'string',para.reg(med).bet);
set(bouton.alpha,'string',para.reg(med).alpha);
set(bouton.rho  ,'string',para.reg(med).rho);

set(bouton.Q    ,'string',para.reg(med).qd);
set(bouton.lstatts  ,'value',para.reg(med).tipoatts);
if para.rafraichi==0
    cmd_att;
    cmd_lambda_mu;
end
if med==1
    set(info.geo    ,'string',['El medio ',num2str(med),' se define como']);
    set(bouton.geo  ,'string',strgeo1,'value',para.geo(med));
    set(info.cont   ,'visible','off');
    set(bouton.cont ,'visible','off');
    
    set(info.DWN,'visible','off');
    set(bouton.nsubmed,'visible','off');
%     set(info.DWNkmax,'visible','off');
%     set(bouton.DWNkmax,'visible','off');
    set(info.DWNnbptkx,'visible','off');
    set(bouton.DWNnbptkx,'visible','off');
    set(info.infoDWN,'visible','off');
    set(info.estrDWN,'visible','off');
    set(info.submed,'visible','off');
    set(bouton.submed,'visible','off');
    set(info.subh,'visible','off');
    set(bouton.subh,'visible','off');
        
    if para.geo(med) == 1
        onoff='off';
        set(bouton.recpos,'value',2);
        para.recpos=2;
    elseif para.geo(med) == 2
        onoff='on';
        set(bouton.xa	,'string',para.cont(med,1).xa);
        set(bouton.a	,'string',para.cont(med,1).a);
        set(bouton.th	,'string',para.cont(med,1).th);
    elseif para.geo(med) == 3
        onoff='off';
        set(bouton.xa	,'string',para.cont(med,1).xa);
        set(bouton.a	,'string',para.cont(med,1).a);
        set(bouton.th	,'string',para.cont(med,1).th);
        
        set(info.DWN,'visible','on');
        set(bouton.nsubmed,'visible','on');
%         set(info.DWNkmax,'visible','on');
%         set(bouton.DWNkmax,'visible','on');
        set(info.DWNnbptkx,'visible','on');
        set(bouton.DWNnbptkx,'visible','on');
        set(info.infoDWN,'visible','on');
        set(info.estrDWN,'visible','on');
        set(info.submed,'visible','on');
        set(bouton.submed,'visible','on');
        set(info.subh,'visible','on');
        set(bouton.subh,'visible','on');
        if para.rafraichi==0
            cmd_nsubmed;
        end
    end
    
    set(info.parpos,'visible',onoff);
    set(info.a      ,'string','longitud(x)','visible',onoff);
    set(bouton.a   	,'visible',onoff);
    set(info.xa     ,'visible',onoff);
    set(bouton.xa   ,'visible',onoff);
    set(info.th     ,'visible',onoff);
    set(bouton.th   ,'visible',onoff);
    
    set(info.partes,'visible','off');
    set(info.contgeo,'visible','off');
    set(bouton.contgeo,'visible','off');
    set(info.za     ,'visible','off');
    set(bouton.za   ,'visible','off');
    set(info.base   ,'visible','off');
    set(bouton.base ,'visible','off');
    set(info.haut  	,'visible','off');
    set(bouton.haut	,'visible','off');
    
    set(info.rug   	,'visible',onoff);
    set(info.ruggeo ,'visible',onoff);
    set(bouton.ruggeo,'visible',onoff);
    
    if para.geo(med) == 2
        set(bouton.ruggeo,'value' ,para.cont(med,1).ruggeo);
        if para.cont(med,1).ruggeo==1
            onoff='off';
        else
            onoff='on';
            icont = get(bouton.cont,'value');
            if icont==2
                para.cont(med,1).rba    = para.cont(med,2).rba;
                para.cont(med,1).rh     = para.cont(med,2).rh;
            end
            set(bouton.rugbase   ,'string',para.cont(med,1).rba);
            set(bouton.rughaut   ,'string',para.cont(med,1).rh);
        end
    end
    set(info.rugbase    ,'visible',onoff);
    set(bouton.rugbase  ,'visible',onoff);
    set(info.rughaut  	,'visible',onoff);
    set(bouton.rughaut	,'visible',onoff);
else
    set(info.geo,'string',['El medio ',num2str(med),' se define por']);
    set(bouton.geo  ,'string',strgeo2,'value',para.geo(med));
    icont = get(bouton.cont,'value');
    if para.geo(med) == 1
        onoff='on';
        set(info.a	,'string','1/2 ancho /x (a)');
        set(info.haut   ,'string','altura del contorno');
        if para.cont(med,icont).geom<6
            onoffb='off';
        else
            onoffb='on';
        end
    elseif para.geo(med) >= 2
        onoff='off';
        set(info.a	,'string','largo(x)');
        set(info.haut   ,'string','espesor de la placa');
    end
    set(info.xa     ,'visible','on');
    set(bouton.xa   ,'visible','on');
    set(info.za     ,'visible','on');
    set(bouton.za   ,'visible','on');
    set(info.a      ,'visible','on');
    set(bouton.a   	,'visible','on');
    
    set(info.partes,'visible','on');
    
    set(info.cont   ,'visible','on');
    set(bouton.cont ,'visible','on');
    
    set(info.contgeo,'visible',onoff);
    set(bouton.contgeo,'visible',onoff);
    set(info.base   ,'visible',onoffb);
    set(bouton.base ,'visible',onoffb);
    set(info.haut   ,'visible','on');
    set(bouton.haut ,'visible','on');
    
    set(info.rug   	,'visible','on');
    set(info.ruggeo ,'visible','on');
    set(bouton.ruggeo,'visible','on');
    
    if para.cont(med,icont).ruggeo==1
        onoff='off';
    else
        onoff='on';
    end
    set(info.rugbase    ,'visible',onoff);
    set(bouton.rugbase  ,'visible',onoff);
    set(info.rughaut  	,'visible',onoff);
    set(bouton.rughaut	,'visible',onoff);
    if para.rafraichi==0
        cmd_cont;
    end
end


