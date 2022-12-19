function cmd_dib0_inc(b_dib,para,bouton,utc,uw,stc,sw)

% if isfield(b_dib(1),'name')
%     if isfield(para,'redraw')
%         tmp=para.redraw;
%         load(b_dib(1).name,'para');
%         para.redraw=tmp;
%     else
%         load(b_dib(1).name,'para');
%     end
% end

%recuperacion de los parametros
ifig = 1;
iinc = 1;
if ishandle(b_dib(1).inc0)
    if (b_dib(1).inc0 ~= 0)
    iinc = get(b_dib(1).inc0,'value');
    end
end
if isfield(b_dib(1)','hold0')
holdtr  = get(b_dib(1).hold0,'value');
else
    if para.redraw==0
        holdtr=2;
    else
        holdtr=1;
    end
end
for i=1:length(b_dib)
    h201(i)   = b_dib(i).h201;
end
prev_common_dessin;
dibujo_iinc0;

