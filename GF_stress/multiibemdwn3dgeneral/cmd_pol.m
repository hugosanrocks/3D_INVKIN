if para.pol==1 && para.dim==1
    set(info.alpha,'visible','off');
    set(bouton.alpha,'visible','off');
    set(info.lambda,'visible','off');
    set(bouton.lambda,'visible','off');
else
    set(info.alpha,'visible','on');
    set(bouton.alpha,'visible','on');
    set(info.lambda,'visible','on');
    set(bouton.lambda,'visible','on');
end

if para.rafraichi==0
    cmd_fuente;
end
if para.dim==1 && (para.geo(1)~=3)
    if para.pol==1
        onoff1='on';
        onoff2='off';
        onoff3='on';
        onoff4='off';
    else
        onoff1='off';
        onoff2='on';
        onoff3='on';
        onoff4='on';
    end
    onoff5='off';
elseif para.dim==1 && para.geo(1)==3
    if para.pol==1
        onoff1='on';
        onoff2='off';
    else
        onoff1='off';
        onoff2='on';
    end
    onoff3='off';
    onoff4='off';
    onoff5='off';
else
    onoff1='on';
    onoff2='on';
    onoff3='off';
    onoff4='off';
    onoff5='on';
end

set(boutonsal.Ux ,'visible',onoff2);
set(boutonsal.Uz ,'visible',onoff2);
% set(boutonsal.Vx ,'visible',onoff2);
% set(boutonsal.Vz ,'visible',onoff2);
set(boutonsal.sxx,'visible',onoff2);
set(boutonsal.szz,'visible',onoff2);
set(boutonsal.sxz,'visible',onoff2);
% set(boutonsal.P  ,'visible',onoff2);
% set(boutonsal.exx,'visible',onoff2);
% set(boutonsal.ezz,'visible',onoff2);
% set(boutonsal.exz,'visible',onoff2);

% set(boutonsal.Ut,'visible',onoff2);
set(boutonsal.UPh,'visible',onoff4);
set(boutonsal.USh,'visible',onoff3);
set(boutonsal.UIh,'visible',onoff3);

set(boutonsal.UPt,'visible',onoff4);
set(boutonsal.USt,'visible',onoff4);


set(boutonsal.Uy ,'visible',onoff1);
% set(boutonsal.Vy ,'visible',onoff1);
set(boutonsal.sxy,'visible',onoff1);
set(boutonsal.syz,'visible',onoff1);

% set(boutonsal.eyy,'visible',onoff1);

set(boutonsal.syy,'visible',onoff5);

% set(boutonsal.exy,'visible',onoff);
% set(boutonsal.eyz,'visible',onoff);