val = get(bouton.rotDespConfGeom,'value');
%axes(bouton.axe_conf_geo)
switch val
  case 1
    rotate3d off
    pan off
  case 2
    pan off
    rotate3d on
  case 3
    rotate3d off
    pan on
end
clear val