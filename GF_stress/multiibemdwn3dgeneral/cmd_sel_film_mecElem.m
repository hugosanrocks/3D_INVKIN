cmd_takingNames;
c_2 = [0.85 0.85 0.9];%[0.3 0.8 1];
c_3 = [0.8 0.8 0.8];
boutonfmeExist = isfield(bouton,'fil_mec');
if para.film.filmStyle == 4 % Mecanical elements on the boundary
  if length(para.film.strmecelemlist) >= 3
  if (strcmp(para.film.strmecelemlist(end-2),'sxx') && ...
      strcmp(para.film.strmecelemlist(end-1),'szz') && ...
      strcmp(para.film.strmecelemlist(end-0),'sxz'))
    para.film.strmecelemlist = [para.film.strmecelemlist 'stt' 'srt'];
  end
  end
  if (boutonfmeExist)
    delete(bouton.fil_mec);
    delete(bouton.fil_BouElemRang);
    delete(bouton.tx_BouElemRang);
    bouton = rmfield(bouton,'fil_mec');
    bouton = rmfield(bouton,'fil_BouElemRang');
    bouton = rmfield(bouton,'tx_BouElemRang');
  end
  para.film.filmeMecElem = 1;
  bouton.fil_mec =uicontrol('parent',bouton.uiwr,'Style','popupmenu','BackgroundColor',c_2,'Units','normalized','position',[0.5 0.34 .2 .1],...
    'string',para.film.strmecelemlist,'Callback','para.film.filmeMecElem = get(bouton.fil_mec,''value'');');
  para.film.BoundaryWarpRange = '1:end';
  bouton.tx_BouElemRang = uicontrol('parent',bouton.uiwr,'Style','text','BackgroundColor',c_3,'Units','normalized','position',[0.71 0.37 .12 .07],...
    'string','Boundaries:');
  bouton.fil_BouElemRang = uicontrol('parent',bouton.uiwr,'Style','edit','BackgroundColor',c_2,'Units','normalized','position',[0.83 0.37 .1 .08],...
    'string',para.film.BoundaryWarpRange,'Callback','para.film.BoundaryWarpRange = get(bouton.fil_BouElemRang,''value'');');
else
  if (boutonfmeExist)
    para.film.filmeMecElem = 1;
    para.film.BoundaryWarpRange = '1:end';
    delete(bouton.fil_mec);
    delete(bouton.fil_BouElemRang);
    delete(bouton.tx_BouElemRang);
    bouton = rmfield(bouton,'fil_mec');
    bouton = rmfield(bouton,'fil_BouElemRang');
    bouton = rmfield(bouton,'tx_BouElemRang');
  end
end
clear c_2 c_3 boutonfmeExist