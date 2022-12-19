% que el indice de tipo de geometría no exceda el tamaño de la lista
if para.dim >= 3
  diNew = get(bouton.dim,'value');
  if (diNew < 3)
    set(bouton.geo,'value',1);
  end
end
drawnow;
para.geo(get(bouton.med,'value'))=get(bouton.geo,'value');
para.chggeo=1;
cmd_med