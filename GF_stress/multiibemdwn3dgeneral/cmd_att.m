tatt=para.reg(get(bouton.med,'value')).tipoatts;
if tatt==1
    set(info.att,'string','Q factor');
    set(bouton.Q,'visible','on');
elseif tatt==2
    set(info.att,'string','KV factor');
    set(bouton.Q,'visible','on');
elseif tatt==0 || tatt==3
    set(info.att,'string','sin atenuacion');
    set(bouton.Q,'visible','off');
end
clear tatt