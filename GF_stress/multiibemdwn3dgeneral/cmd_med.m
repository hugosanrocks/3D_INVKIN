med     =get(bouton.med,'value');

if para.chggeo==1
%   if para.smallscreen
    para.cont(med,1).piece{info.ThisPiece}.fileName = '';
    para.cont(med,1).piece{info.ThisPiece}.geoFileData = [];
%   end
  % Nuevo elemento de geometría. Inicializar variables.
    para.chggeo=0;
    if para.rafraichi==0
        cmd_pol;%chgmt des sorties possibles
    end
    if med==1 && para.geo(1)==3 %&& para.dim <4
      % fondo N-layer DWN
        if para.nsubmed>2
            %il faut supprimer des champs;
            sub0            = para.reg(1).sub(1:2);
            para.reg(1).sub = sub0;
        end
        para.nsubmed   = 2;
        for i=1:para.nsubmed
            para.reg(1).sub(i).rho      = 1;
            para.reg(1).sub(i).alpha    = 2;
            para.reg(1).sub(i).bet      = 1;
            para.reg(1).sub(i).qd       = 1000;
            para.reg(1).sub(i).tipoatts = 1;
            para.reg(1).sub(i).lambda   = para.reg(1).sub(i).rho*(para.reg(1).sub(i).alpha^2-2*para.reg(1).sub(i).bet^2);
            para.reg(1).sub(i).mu       = para.reg(1).sub(i).rho*para.reg(1).sub(i).bet^2;
            para.reg(1).sub(i).h        = 1;
        end
        para.reg(1).sub(i).h        = 0;%ultimo estrato=semi espacio
        strsubmed      = 1:para.nsubmed;
        set(bouton.submed ,'string',strsubmed,'value',1);
        set(bouton.nsubmed,'string',para.nsubmed);
        
    elseif para.geo(med) == 1 && med>1 && para.dim <4
      % inclusión: 2 contornos horizontales
        para.cont(med,1).ba   = 0.25;
        para.cont(med,2).ba   = 0.25;
    elseif para.geo(med) == 2 && med>1 && para.dim <4
      % inclusión: placa ilimitada
        para.cont(med,1).ba   = 0;
        para.cont(med,2).ba   = 0;
        if para.geo(1)==2
            para.cont(med,1).a =para.cont(1,1).a;
            para.cont(med,1).xa=para.cont(1,1).xa;
        else
            para.cont(med,1).a = 6;
            para.cont(med,1).xa=-3;
        end
    elseif (para.geo(med) == 3 || para.geo(med) == 4) && med>1 && para.dim <4
      % inclusión: semiplaca Left o Right
        para.cont(med,1).ba   = 0.2;
        para.cont(med,2).ba   = 0.3;
        if para.geo(1)==2
            para.cont(med,1).a =para.cont(1,1).a/2;
            if  para.geo(med)==3
                para.cont(med,1).xa=0;
            else
                para.cont(med,1).xa=para.cont(1,1).xa;
            end
        else
            para.cont(med,1).a = 6;
            para.cont(med,1).xa=-3;
        end
    elseif  para.dim == 4 && med > 1
%       if para.smallscreen
      % inclusión 3D: desde archivo STL
      para.cont(med,1).piece{info.ThisPiece}.fileName = get(bouton.geoFileSelect,'string');
      if ~strcmp(para.cont(med,1).piece{info.ThisPiece}.fileName,'')
        [~,~] = previewSTL(bouton.gfPreview,para.cont(med,1).piece{info.ThisPiece});
      end
%       end
    end
end

% Actualziar los botones en "Propiedades de los medios"
if med==1 && para.geo(1) == 3
    % fondo N-layer DWN
    % il faut changer les callback des champs rho,beta,alpha,Q
    set(bouton.rho    ,'Callback','para.reg(1).sub(get(bouton.submed,''value'')).rho=str2num(get(bouton.rho,''string''));cmd_lambda_mu_sub;');
    set(bouton.alpha  ,'Callback','para.reg(1).sub(get(bouton.submed,''value'')).alpha=str2num(get(bouton.alpha,''string''));cmd_lambda_mu_sub;');
    set(bouton.bet    ,'Callback','para.reg(1).sub(get(bouton.submed,''value'')).bet=str2num(get(bouton.bet,''string''));cmd_lambda_mu_sub;');
    set(bouton.lstatts,'Callback','para.reg(1).sub(get(bouton.submed,''value'')).tipoatts=get(bouton.lstatts,''value'');cmd_att_sub;');
    set(bouton.Q      ,'Callback','para.reg(1).sub(get(bouton.submed,''value'')).qd =str2num(get(bouton.Q  ,''string''));');
    
    set(bouton.bet  ,'string',para.reg(1).sub(1).bet);
    set(bouton.alpha,'string',para.reg(1).sub(1).alpha);
    set(bouton.rho  ,'string',para.reg(1).sub(1).rho);
    set(bouton.Q    ,'string',para.reg(1).sub(1).qd);
    set(bouton.lstatts  ,'value',para.reg(1).sub(1).tipoatts);
    if para.rafraichi==0
        cmd_att_sub;
        cmd_lambda_mu_sub;
    end
    strsubmed   = 1:para.nsubmed;
    set(bouton.submed,'string',strsubmed,'value',1);
    set(bouton.nsubmed,'string',para.nsubmed);
    if para.rafraichi==0
        cmd_nsubmed;
    end
else
  % fondo con IBEM ó una inclusión la que sea
    % il faut changer les callback des champs rho,beta,alpha,Q si med==1 && geo==3
    set(bouton.rho    ,'Callback','para.reg(get(bouton.med,''value'')).rho=str2num(get(bouton.rho,''string''));cmd_lambda_mu;');
    set(bouton.alpha  ,'Callback','para.reg(get(bouton.med,''value'')).alpha=str2num(get(bouton.alpha,''string''));cmd_lambda_mu;');
    set(bouton.bet    ,'Callback','para.reg(get(bouton.med,''value'')).bet=str2num(get(bouton.bet,''string''));cmd_lambda_mu;');
    set(bouton.lstatts,'Callback','para.reg(get(bouton.med,''value'')).tipoatts=get(bouton.lstatts,''value'');cmd_att;');
    set(bouton.Q      ,'Callback','para.reg(get(bouton.med,''value'')).qd =str2num(get(bouton.Q  ,''string''));');
    
    set(bouton.bet  ,'string',para.reg(med).bet);
    set(bouton.alpha,'string',para.reg(med).alpha);
    set(bouton.rho  ,'string',para.reg(med).rho);
    set(bouton.Q    ,'string',para.reg(med).qd);
    set(bouton.lstatts  ,'value',para.reg(med).tipoatts);
    if para.rafraichi==0
        cmd_att;
        cmd_lambda_mu;
    end
end

% Actualizar los botones en "Geometría del medio"
if med==1
  % Del fondo
    set(info.geo    ,'string',['El medio ',num2str(med),' se define como']);
    set(bouton.geo  ,'string',bouton.strgeo1,'value',para.geo(med));
    set(info.cont   ,'visible','off');
    set(bouton.cont ,'visible','off');
    
    set(info.DWN,'visible','off');
    set(bouton.nsubmed,'visible','off');
    set(info.DWNxl,'visible','off');
    set(bouton.DWNxl,'visible','off');
    set(info.DWNnbptkx,'visible','off');
    set(bouton.DWNnbptkx,'visible','off');
%     set(info.infoDWN,'visible','off');
    set(info.estrDWN,'visible','off');
    set(info.submed,'visible','off');
    set(bouton.submed,'visible','off');
    set(info.subh,'visible','off');
    set(bouton.subh,'visible','off');
    
    if para.geo(1) == 1
        onoff='off';
        if para.recpos==1
            set(bouton.recpos,'value',2);
            para.recpos=2;
        end
    elseif para.geo(1) == 2
        onoff='on';
        set(bouton.xa	,'string',para.cont(1,1).xa);
        set(bouton.a	,'string',para.cont(1,1).a);
        set(bouton.th	,'string',para.cont(1,1).th);
    elseif para.geo(1) == 3
        onoff='off';
        set(info.DWN,'visible','on');
        set(bouton.nsubmed,'visible','on');
        set(info.DWNxl,'visible','on');
        set(bouton.DWNxl,'visible','on');
        set(info.DWNnbptkx,'visible','on');
        set(bouton.DWNnbptkx,'visible','on');
%         set(info.infoDWN,'visible','on');
        set(info.estrDWN,'visible','on');
        set(info.submed,'visible','on');
        set(bouton.submed,'visible','on');
        set(info.subh,'visible','on');
        set(bouton.subh,'visible','on');
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
%    if para.smallscreen
    set(bouton.geoFileSelect,'visible','off');
    set(bouton.geoFilePanel,'visible','off');
    set(allchild(bouton.gfPreview),'visible','off')
    if para.geo(1) == 2
        set(bouton.ruggeo,'value' ,para.cont(1,1).ruggeo);
        if para.cont(1,1).ruggeo==1
            onoff='off';
        else
            onoff='on';
            icont = get(bouton.cont,'value');
            if icont==2
                para.cont(1,1).rba    = para.cont(1,2).rba;
                para.cont(1,1).rh     = para.cont(1,2).rh;
            end
            set(bouton.rugbase   ,'string',para.cont(1,1).rba);
            set(bouton.rughaut   ,'string',para.cont(1,1).rh);
        end
    end
%    end
    set(info.rugbase    ,'visible',onoff);
    set(bouton.rugbase  ,'visible',onoff);
    set(info.rughaut  	,'visible',onoff);
    set(bouton.rughaut	,'visible',onoff);
else
  % De una inclusión
    set(info.geo,'string',['El medio ',num2str(med),' se define por']);
    if (para.dim <= 2) % en 2D ó 2.5D
    set(bouton.geo  ,'string',bouton.strgeo2,'value',para.geo(med));
    elseif (para.dim == 3) % 3D axisimétrico
    set(bouton.geo  ,'string',bouton.strgeo3,'value',para.geo(med));
    elseif (para.dim == 4) % 3D geometría arbitraria desde archivos
    set(bouton.geo  ,'string',bouton.strgeo4,'value',para.geo(med));
    end
    icont = get(bouton.cont,'value');
    onoffc='on'; onoffd = 'off';
   if  para.dim <4
    if para.geo(med) == 1
        %contours fermés constitués de 2 sous contours
        onoff='on';
        set(info.a	,'string','1/2 ancho /x (a)');
        set(info.haut   ,'string','altura del contorno');
        if para.cont(med,icont).geom<6
            onoffb='off';
        else
            onoffb='on';
        end
    elseif para.geo(med) == 2
        %plaque infinie
        onoff='off';
        set(info.a	,'string','largo(x)');
        set(info.haut   ,'string','espesor de la placa');
        onoffb='off';
    elseif para.geo(med) == 3 || para.geo(med) == 4
        %semi plaque infinie a droite ou a gauche
        onoff='on';
        set(info.a	,'string','largo(x)');
        set(info.haut   ,'string','espesor de la placa');
        onoffb='on';
    end
   elseif para.dim == 4
        % partir d'un fichier STL
        onoff = 'off';
        onoffb = 'off';
        onoffc = 'off';
        onoffd = 'on';
    end
    set(info.parpos  ,'visible',onoffc);
    set(info.xa     ,'visible',onoffc);
    set(bouton.xa   ,'visible',onoffc);
    set(info.za     ,'visible',onoffc);
    set(bouton.za   ,'visible',onoffc);
    set(info.a      ,'visible',onoffc);
    set(bouton.a   	,'visible',onoffc);
    set(info.th     ,'visible',onoffc);
    set(bouton.th   ,'visible',onoffc);
    set(info.partes,'visible',onoffc);
    
    set(info.cont   ,'visible',onoffc);
    set(bouton.cont ,'visible',onoffc);
    
    set(info.contgeo,'visible',onoff);
    set(bouton.contgeo,'visible',onoff);
    set(info.base   ,'visible',onoffb);
    set(bouton.base ,'visible',onoffb);
    set(info.haut   ,'visible',onoffc);
    set(bouton.haut ,'visible',onoffc);
    
    set(info.rug   	,'visible',onoffc);
    set(info.ruggeo ,'visible',onoffc);
    set(bouton.ruggeo,'visible',onoffc);
    set(bouton.geoFileSelect,'visible',onoffd);
    set(bouton.geoFilePanel,'visible',onoffd);
    set(allchild(bouton.gfPreview),'visible',onoffd);
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
    
    set(info.DWN,'visible','off');
    set(bouton.nsubmed,'visible','off');
    set(info.DWNxl,'visible','off');
    set(bouton.DWNxl,'visible','off');
    set(info.DWNnbptkx,'visible','off');
    set(bouton.DWNnbptkx,'visible','off');
%     set(info.infoDWN,'visible','off');
    set(info.estrDWN,'visible','off');
    set(info.submed,'visible','off');
    set(bouton.submed,'visible','off');
    set(info.subh,'visible','off');
    set(bouton.subh,'visible','off');
end
if para.rafraichi==0
    cmd_fuente;%pour update des ondes planes possibles si med==1 && para.geo(1)==3
end
% drawnow

% clear alpha bet icont lambda mu onoff onoffb onoffc
% clear onoffd rho strtipo_onda tatt