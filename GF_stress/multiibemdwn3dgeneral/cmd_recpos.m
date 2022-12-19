if para.dim >= 3
  if para.recpos == 1; 
    para.recpos = 2;
  end
end

set(bouton.recpos   ,'value'    ,para.recpos);
% if para.smallscreen
set(bouton.resatboundary ,'value', para.rec.resatboundary);
set(bouton.resatboundaryDecimate, 'string', para.rec.resatboundaryDecimate);
set(bouton.resatboundaryScale, 'string', para.rec.resatboundaryScale);
set(bouton.resatboundaryScaleMediaRange, 'string',para.rec.resatboundaryScaleMediaRange);
% end
set(bouton.xri      ,'string'   ,para.rec.xri);
set(bouton.dxr      ,'string'   ,para.rec.dxr);
set(bouton.nrecx    ,'string'   ,para.rec.nrecx);
set(bouton.yri      ,'string'   ,para.rec.yri);
set(bouton.dyr      ,'string'   ,para.rec.dyr);
set(bouton.nrecy    ,'string'   ,para.rec.nrecy);
set(bouton.zri      ,'string'   ,para.rec.zri);
set(bouton.dzr      ,'string'   ,para.rec.dzr);
set(bouton.nrecz    ,'string'   ,para.rec.nrecz);

if para.recpos==1
    %superficie
    set(info.recdr  ,'visible','on');
    set(bouton.dxr  ,'visible','on');
    set(bouton.xri  ,'callback','para.rec.xri=str2num(get(bouton.xri,''string''));');
    
    if para.dim==1
        para.rec.nrecy=1;
        set(info.recy   ,'visible','off');
        set(bouton.yri  ,'visible','off');
        set(bouton.dyr  ,'visible','off');
        set(bouton.nrecy,'visible','off');
    else
        set(info.recy   ,'visible','on');
        set(bouton.yri  ,'visible','on');%,'callback','para.rec.yri=str2num(get(bouton.yri,''string''));');
        set(bouton.dyr  ,'visible','on');
        set(bouton.nrecy,'visible','on');
    end
    
    para.rec.nrecz=1;
    set(info.recz   ,'visible','off');
    set(bouton.zri  ,'visible','off');
    set(bouton.dzr  ,'visible','off');
    set(bouton.nrecz,'visible','off');

    set(info.irec,'visible','off');
    set(bouton.irec,'visible','off');
elseif para.recpos==2
    %malla
    set(info.recdr  ,'visible','on');
    set(bouton.dxr  ,'visible','on');
    if para.dim<3
        set(info.recy   ,'visible','off');
        set(bouton.yri  ,'visible','off');
        set(bouton.dyr  ,'visible','off');
        set(bouton.nrecy,'visible','off');
    else
        set(info.recy   ,'visible','on');
        set(bouton.yri  ,'visible','on');
        set(bouton.dyr  ,'visible','on');
        set(bouton.nrecy,'visible','on');
    end
    set(info.recz   ,'visible','on');
    set(bouton.zri  ,'visible','on');
    set(bouton.dzr  ,'visible','on');
    set(bouton.nrecz,'visible','on');
    
    set(info.irec,'visible','off');
    set(bouton.irec,'visible','off');
    set(bouton.xri  ,'callback','para.rec.xri=str2num(get(bouton.xri,''string''));');
    set(bouton.yri  ,'callback','para.rec.yri=str2num(get(bouton.yri,''string''));');
    set(bouton.zri  ,'callback','para.rec.zri=str2num(get(bouton.zri,''string''));');
elseif para.recpos==3
    %posicion libre
    if para.chgrec==1
        para.chgrec=0;
        para.rec.xr=zeros(para.rec.nrecx,1);
        para.rec.zr=zeros(para.rec.nrecx,1);
        para.rec.yr=zeros(para.rec.nrecx,1);
    end
    if para.chgnrecx==1
        para.chgnrecx=0;
        strrec	= 1:para.rec.nrecx;
        set(bouton.irec,'string',strrec,'value',1);
        if para.rec.nrecx>=length(para.rec.xr)
            %il faut initialiser des champs;
            i=length(para.rec.xr)+1:para.rec.nrecx;
            para.rec.xr(i)=0;
            para.rec.zr(i)=0;
            para.rec.yr(i)=0;
        else
            %il faut supprimer des champs;
            xr          = para.rec.xr(1:para.rec.nrecx);
            yr          = para.rec.yr(1:para.rec.nrecx);
            zr          = para.rec.zr(1:para.rec.nrecx);
            para.rec.xr = xr;
            para.rec.yr = yr;
            para.rec.zr = zr;
        end
    end
    
    set(info.recdr  ,'visible','off');
    set(bouton.dxr  ,'visible','off');
    set(bouton.dyr  ,'visible','off');
    set(bouton.nrecy,'visible','off');
    
    set(info.recz   ,'visible','on');
    set(bouton.zri  ,'visible','on');
    set(bouton.dzr  ,'visible','off');
    set(bouton.nrecz,'visible','off');
    
    set(info.irec,'visible','on');
    set(bouton.irec,'visible','on');
    irec = get(bouton.irec,'value');
    set(bouton.xri  ,'string',para.rec.xr(irec),'callback',['para.rec.xr(',num2str(irec),')=str2num(get(bouton.xri,''string''));']);
    set(bouton.zri  ,'string',para.rec.zr(irec),'callback',['para.rec.zr(',num2str(irec),')=str2num(get(bouton.zri,''string''));']);
    if para.dim<3
        set(info.recy   ,'visible','off');
        set(bouton.yri  ,'visible','off');
    else
        set(info.recy   ,'visible','on');
        set(bouton.yri  ,'visible','on','string',para.rec.yr(irec),'callback',['para.rec.yr(',num2str(irec),')=str2num(get(bouton.yri,''string''));']);
    end
end