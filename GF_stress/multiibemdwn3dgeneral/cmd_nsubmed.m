nsubmednew=str2double(get(bouton.nsubmed,'string'));

if nsubmednew>para.nsubmed
    %il faut initialiser des champs;
    for i=para.nsubmed+1:nsubmednew
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
elseif nsubmednew<para.nsubmed
    %il faut supprimer des champs;
    if nsubmednew==0
        nsubmednew=1;
        warndlg('No seas ...');
    end
    sub0            = para.reg(1).sub(1:nsubmednew);
    para.reg(1).sub = sub0;
    para.reg(1).sub(nsubmednew).h        = 0;%ultimo estrato=semi espacio
    
end

para.nsubmed   = nsubmednew;
strsubmed      = 1:para.nsubmed;
set(bouton.submed,'string',strsubmed,'value',1);
cmd_submed;