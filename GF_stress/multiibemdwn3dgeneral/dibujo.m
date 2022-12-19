function para = dibujo(para,bouton,RESULT)
utc = RESULT.utc;
uw = RESULT.uw;
stc = RESULT.stc;
sw = RESULT.sw;
% cont1 = RESULT.cont1;

if ~isstruct(para.b_dib)
    clear para.b_dib
end

%recuperacion y inicialisacion de algunos parametros 
ninc    = para.ninc;

prev_common_dessin;
para.b_dib(1).dessinsp   = get(bouton.dessinsp,'value');
para.b_dib(1).normalise  = get(bouton.normalise,'value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figura con senales en tiempo y geometria %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h201 = zeros(ns0,1);
for ifig=1:ns0 %numero de campos
    %pour chacune des sorties genérales Ut, UP, US, UI, ou S (hors composantes)
    %(a etendre a sigma, v, def), on cree une figure
    h201(ifig)=figure(200+ifig);
    clf
    set(h201(ifig),'numbertitle','off');
    
    para.b_dib(ifig).hold0=  uicontrol('Style','popupmenu','BackgroundColor',[0.3 0.8 1],'Units','normalized','position',[0.01 0.32 .07 .05],...
            'string',{'on','off'}); 
     if para.ninc > 1
        %bouton pour superposer ou non les signaux
        uicontrol('Style','text'     ,'BackgroundColor',[1 .1 .1]  ,'Units','normalized','position',[0.01 0.37 .07 .05],'string','hold');
        
   
        %bouton pour choisir 1 signal parmis tous
        strinc =[];
        for i=1:ninc
            tmpstr=num2str(i);
            while length(tmpstr)<floor(log10(ninc)+1)
                tmpstr=['0',tmpstr];
            end
            tmpstr = ['Inc_',tmpstr];
            strinc =[strinc;tmpstr];
        end
                            uicontrol('Style','text'     ,'BackgroundColor',[1 .1 .1]  ,'Units','normalized','position',[0.01 0.17 .07 .05],'string','Inc. #');
        para.b_dib(ifig).inc0  =  uicontrol('Style','popupmenu','BackgroundColor',[0.3 0.8 1],'Units','normalized','position',[0.01 0.12 .07 .05],...
            'string',strinc,'Callback',['ifig=',num2str(ifig),';cmd_dib0_inc(para.b_dib,para,bouton,RESULT.utc,RESULT.uw,RESULT.stc,RESULT.sw);']);
    else
%         set(para.b_dib(ifig).hold0,'Visible','off');
        para.b_dib(ifig).inc0  = 0;
    end
    %button to chose the wiggle format
        para.b_dib(ifig).wiggle = uicontrol('Style','popupmenu',...
            'BackgroundColor',[0.3 0.8 1],...
            'Units','normalized','position',[0.7 0.95 .2 .05],...
            'string',{'wiggle on','wiggle off'},'value',2,...
            'Callback','cmd_dib0_inc(para.b_dib,para,bouton,RESULT.utc,RESULT.uw,RESULT.stc,RESULT.sw);');
    %panel para que la figura no se obstruya con los botones
        uipan = uipanel('parent',gcf,'Position',[0 0 1 .95]);
        axes('parent',uipan)
        
    %%%%%%%%%%%%%
    % geometria %
    %%%%%%%%%%%%%
    nxyz=fieldV(ifig).nc;
    para.b_dib(ifig).hg=subplot(1,(nxyz+2),1);
    
    hold on
    xlabel('x')
    ylabel('z')
    
    nmed1  = length(RESULT.cont1);
    for i = 1:nmed1
        plot(RESULT.cont1(i).vec.xc,RESULT.cont1(i).vec.zc,'r')
        if para.dim>=3
            plot(-RESULT.cont1(i).vec.xc,RESULT.cont1(i).vec.zc,'r')
        end
    end
    
    %receptores
    plot(xr,zr,'.b');
    
    set(gca,'ydir','reverse');
    if para.rec.nrecz>1 && para.rec.nrecx==1
    else
        view(-90,90)
    end
    tmp=get(para.b_dib(ifig).hg,'title');
    set(tmp,'string','geometria');
    
    if para.spct==0
        %%%%%%%%%%%%%%%%%%%%%
        % senales en tiempo %
        %%%%%%%%%%%%%%%%%%%%%
        pos=get(para.b_dib(ifig).hg,'position');
        x0hs=1.2*(pos(1)+pos(3));%X0 des autres graphs
        
        DX = (1-x0hs)/(fieldV(ifig).nc);
        for ic=1:fieldV(ifig).nc
            para.b_dib(ifig).ht(ic)=subplot('Position',[x0hs+(ic-1)*DX pos(2) DX/1.3 pos(4)]) ;hold on;
            xlabel(para.b_dib(ifig).ht(ic),'t en (s)')
            tmp  = get(para.b_dib(ifig).ht(ic),'title');
            tmpc0= fieldV(ifig).namec{ic};
            set(tmp,'string',[fieldV(ifig).name,'_{',tmpc0(2:end),'}']);
        end

        hold on;
        if para.redraw==0
            holdtr=2;
        else
            holdtr=1;
        end
        set(h201(ifig),'Toolbar','figure');
        iinc = 1;
        if  isfield(para.b_dib(ifig),'h201')
            if ishandle(para.b_dib(ifig).h201)
                if ishandle(para.b_dib(ifig).inc0)
                    if isfield (para.b_dib(ifig).inc0,'value')
                        iinc = get(para.b_dib(ifig).inc0,'value');
                    end
                end
            end
        end
        dibujo_iinc0;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% senal en tiempo y frecuencia de un receptor %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ifig=1:ns0
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
            
    h101(ifig)=figure(100+ifig);
    set(h101(ifig),'numbertitle','off');
        %bouton pour superposer ou non les signaux
%         ;       	uicontrol('Style','text'     ,'BackgroundColor',[1 .1 .1]  ,'Units','normalized','position',[0.93 0.37 .07 .05],'string','hold');
        para.b_dib(ifig).hold= uicontrol('Style','popupmenu','BackgroundColor',[0.3 0.8 1],'Units','normalized','position',[0.1 0.95 .2 .05],...
            'string',{'hold on','hold off'},'value',1); 
    if (para.rec.nrec==1)
        set(para.b_dib(ifig).hold,'Visible','off');
    end
    

    if ninc>1
        % selecion de la fuente OP/FP# i
%         ;            uicontrol('Style','text'     ,'BackgroundColor',[1 .1 .1]  ,'Units','normalized','position',[0.93 0.17 .07 .05],'string','Inc. #');
        para.b_dib(ifig).inc  = uicontrol('Style','popupmenu','BackgroundColor',[0.3 0.8 1],'Units','normalized','position',[0.3 0.95 .2 .05],...
            'string',strinc,'Callback',['whos bouton;ifig=',num2str(ifig),';cmd_dib_inc_rec(para,bouton,RESULT.utc,RESULT.uw,RESULT.stc,RESULT.sw,ifig);'],'value',iinc);
    else
        para.b_dib(ifig).inc  = 0;
    end
    
    if nrec>1
        strnrec =[];
        for i=1:nrec
            tmpstr=num2str(i);
            while length(tmpstr)<floor(log10(nrec)+1)
                tmpstr=['0',tmpstr];
            end
            tmpstr = ['Res_',tmpstr];
            strnrec =[strnrec;tmpstr];
        end
        % selecion del receptor # i
%         ;            uicontrol('Style','text'     ,'BackgroundColor',[1 .1 .1]  ,'Units','normalized','position',[0.93 0.06 .07 .05],'string','Rec. #');
        
        para.b_dib(ifig).nrec = uicontrol('Style','popupmenu','BackgroundColor',[0.3 0.8 1],'Units','normalized','position',[0.5 0.95 .2 .05],...
            'string',strnrec,'value',irec,'Callback',['ifig=',num2str(ifig),';cmd_dib_inc_rec(para,bouton,RESULT.utc,RESULT.uw,RESULT.stc,RESULT.sw,ifig);']);
        
    else
        para.b_dib(ifig).nrec  = 0;
    end
    
    if para.redraw==0
        holdtr=2;
    else
        holdtr=1;
    end
    set(h101(ifig),'Toolbar','figure');
    %panel para que la figura no se obstruya con los botones
    uipan = uipanel('parent',gcf,'Position',[0 0 1 .95]);
    axes('parent',uipan);
    dibujo_iinc_irec;
end

for ifig=1:ns0
    para.b_dib(ifig).h201=h201(ifig);
    para.b_dib(ifig).h101=h101(ifig);
end


%avance del color
coul    = get(bouton.couleur,'value');
set(bouton.couleur,'value',mod(coul,length(tmpc))+1);
end