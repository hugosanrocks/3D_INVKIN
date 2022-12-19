if para.dim>1
%     set(info.polst,'visible','off');
    set(bouton.pol,'visible','off');
    set(info.phi  ,'visible','on');
    set(bouton.phi,'visible','on');
    set(info.ys   ,'visible','on');
    set(bouton.ys ,'visible','on');
else
%     set(info.polst,'visible','on');
    set(bouton.pol,'visible','on');
    set(info.phi  ,'visible','off');
    set(bouton.phi,'visible','off');
    set(info.ys   ,'visible','off');
    set(bouton.ys ,'visible','off');
end
if para.rafraichi==0
    cmd_pol;
    cmd_recpos;
end