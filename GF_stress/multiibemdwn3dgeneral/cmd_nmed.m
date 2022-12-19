nmednew=str2double(get(bouton.nmed,'string'));

if nmednew>para.nmed
    %il faut initialiser des champs;
    for i=para.nmed+1:nmednew
        para.reg(i).rho     = 1;
        para.reg(i).alpha   = 1;
        para.reg(i).bet     = 1;
        para.reg(i).qd      = 1000;
        para.reg(i).tipoatts= 1;
        para.reg(i).lambda	= para.reg(i).rho*(para.reg(i).alpha^2-2*para.reg(i).bet^2);
        para.reg(i).mu      = para.reg(i).rho*para.reg(i).bet^2;
        para.geo(i)         = 1;
        if nmednew==2
            para.cont(2,1).a = 1;
            para.cont(2,1).th= 0;
            para.cont(2,1).xa=-1;
            para.cont(2,1).za= 0;
            for j=1:2
                para.cont(2,j).geom = 7;
                para.cont(2,j).ba   = 0.6;
            end
            para.cont(2,1).h	=-0.2;
            para.cont(2,2).h	= 0.5;
        else
            para.cont(i,1).a = .25;
            para.cont(i,1).th= 0;
            for j=1:2
                para.cont(i,j).geom = 5;
                para.cont(i,j).h	=(-(j==1)+(j==2))*para.cont(i,1).a;
                para.cont(i,j).ba   = 0.6;
            end
            para.cont(i,1).xa=-0.5+i*.5;
            para.cont(i,1).za=-0.5+i*.5;
        end
        for j=1:2
            para.cont(i,j).ruggeo   = 1;
            para.cont(i,j).rba      = 0.1;
            para.cont(i,j).rh       = 0.1;
        end
        para.cont(i,1).piece{info.ThisPiece}.fileName = '';
        para.cont(i,1).piece{info.ThisPiece}.geoFileData = [];
        para.cont(i,1).piece{info.ThisPiece}.kind = 0; % 0: No determinado por este parametro
    end
elseif nmednew<para.nmed
     if nmednew==0
        nmednew=1;
        warndlg('No seas ...'); % lol !
    end
    %il faut supprimer des champs;
    cont0       = para.cont(1:nmednew,:);
    reg0        = para.reg(1:nmednew);
    geo0        = para.geo(1:nmednew);
    para        = rmfield(para,'reg');
    para        = rmfield(para,'cont');
    para        = rmfield(para,'geo');
    para.cont   = cont0;
    para.reg    = reg0;
    para.geo    = geo0;
end

para.nmed   = nmednew;
strmed      = 1:para.nmed;
set(bouton.med,'string',strmed,'value',1);

clear i j nmednew strmed strtipo_onda