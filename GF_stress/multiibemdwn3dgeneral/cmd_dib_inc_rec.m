function cmd_dib_inc_rec(para,bouton,utc,uw,stc,sw,ifig)

% if isfield(para.b_dib(1),'name')
%     if isfield(para,'redraw')
%         tmp=para.redraw;
%         load(para.b_dib(1).name,'para');
%         para.redraw=tmp;
%     else
%         load(para.b_dib(1).name,'para');
%     end
% end
irec = 1;
iinc = 1;
if  ishghandle(100+ifig) && isstruct(para.b_dib)
    %si il existe deja une figure avec un recepteur/incidence en cours
    %conserver le numero courant
    if para.rec.nrec>1
        if isfield(para.b_dib(ifig),'nrec')
            if ishandle(para.b_dib(ifig).nrec)
                irec = get(para.b_dib(ifig).nrec,'value');
                if irec==0;irec=1;end
            end
        end
    end
    if para.ninc>1
        if isfield(para.b_dib(ifig),'inc')
            if ishandle(para.b_dib(ifig).inc)
                iinc = get(para.b_dib(ifig).inc,'value');
                if iinc==0;iinc=1;end
            end
        end
    end
    if irec>para.rec.nrec
        irec=1;
    end
    if iinc>para.ninc
        iinc=1;
    end
else
    if isfield(para.b_dib(ifig),'nrec')
        para.b_dib(ifig).nrec=0;
    end
    if isfield(para.b_dib(ifig),'inc')
        para.b_dib(ifig).inc=0;
    end
end
holdtr  = get(para.b_dib(ifig).hold,'value');

for i=1:length(para.b_dib)
    h101(i)   = para.b_dib(i).h101;
end
prev_common_dessin;
dibujo_iinc_irec;