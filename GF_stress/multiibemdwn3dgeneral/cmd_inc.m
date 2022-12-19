i=get(bouton.inc,'value');
set(bouton.gam      ,'string',para.gam(i));
set(bouton.xs       ,'string',para.xs( i));
set(bouton.zs       ,'string',para.zs( i));
if para.fuente==1
    set(bouton.tipo_onda,'value',para.tipo_onda( i));
end
if para.dim>=3
    set(bouton.ys  	,'string',para.ys( i));
    set(bouton.phi	,'string',para.phi(i));
end

% if para.rafraichi==0
%     dibujo_conf_geo(para,bouton.axe_conf_geo,bouton);
% end